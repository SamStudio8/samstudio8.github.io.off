---
layout: post
title: "Aligned Annihilation II: Dumpster Diving"
excerpt: "I tried to extract a single integer from a core dump and instead fell in to an abyss and learned how to be a computer."
---

<p class="message"><b>Warning</b></br>This is going to get a bit technical and probably quite boring.<p>


If you haven't already looked at what led me to this awful situation, check out what happened when
I [annihilated all my alignment]({% post_url 2015-05-01-aligned-annihilation %}) data.

I've managed to get `rapsearch` to generate two core dumps for inspection with `gdb` but I am struggling to
extract any useful information beyond what I already know. The `abort()` and `raise()` can clearly be
seen in the backtrace and prior to that, construction of a standard library exception. I presume `rapsearch`
was not compiled with debugging symbols, leaving us with those `??` frames.

{% gist 8c88fdd06cd250872d67 backtrace.gdb %}

Desperately attempting to avoid having to edit and recompile `rapsearch`, I began nosing around the core
in a similar fashion as one would poke a stick around in a dirty pond. At first I naively tried to explore
frame 3, treating `0x31926bcbd6` as "the exception" before realising the address was for a function.
If we translate ("unmangle"[^8]) the symbol we can guess it is responsible
for handling assignment of an exception[^7]:

```
std::__exception_ptr::exception_ptr::operator=(std::__exception_ptr::exception_ptr const&)
```

I'd have to try harder.

I started looking at the registers for the same frame, as the function accepts an exception pointer as a parameter
it should be stored in a register here. It was when I started reading about [x86 Calling Conventions](http://en.wikipedia.org/wiki/X86_calling_conventions#System_V_AMD64_ABI) that I realised I was
probably in over my head here. But I was victim to the sunken cost fallacy now, I had to continue.
I tried to disassemble the frame:

{% gist 8c88fdd06cd250872d67 disassemble_f3.gdb %}

It's [my understanding](http://stackoverflow.com/questions/2535989/what-are-the-calling-conventions-for-unix-linux-system-calls-on-x86-64)
that the `%rdi` register contains the first parameter given to a function and so here I expect to see
some form of pointer to exception. I'll dump the values held in the registers at this frame too.

{% gist 8c88fdd06cd250872d67 registers_f3.gdb %}

Unfortunately `0xb96` doesn't appear to point to anything. I wondered whether it would be
worth trying to follow through the assembly instructions instead.
After determining [which way around](http://en.wikibooks.org/wiki/X86_Assembly/GAS_Syntax) to even read
the instructions[^9], I could become my own arch-nemesis, a computer:

```
[...]
sub    $0x28,%rsp   # Substract 28 from the %rsp in-place
                    # 0x7ffffffbb3e0 − 0x28 = 0x7ffffffbb3b8
                    # %rsp = 0x7ffffffbb3b8
[...]
mov    %rsp,%rdi    # Move %rsp to %rdi
                    # %rdi = 0x7ffffffbb3b8
callq  0x31926559d8 <_ZNSt15__exception_ptr13exception_ptrC1ERKS0_@plt>
```
Ok, so that's a start. `_ZNSt15__exception_ptr13exception_ptrC1ERKS0_` unmangles[^8] to:

```
std::__exception_ptr::exception_ptr::exception_ptr(std::__exception_ptr::exception_ptr const&)
```

A constructor! Expecting a reference to pointer as its first parameter. So what's did we leave `%rdi`
after the subtraction? It looks address-worthy, let's e**x**amine it:

```
(gdb) x 0x7ffffffbb3b8
0x7ffffffbb3b8: 0x23b8db80
(gdb) x 0x23b8db80
0x23b8db80:     0x00463b40
(gdb) x 0x00463b40
0x463b40 <_ZTIN5boost7archive17archive_exceptionE>:     0x0069ee50
```

Seems promising? We're hunting for information on an [`archive_exception`](http://www.boost.org/doc/libs/1_37_0/libs/serialization/doc/exceptions.html)! `_ZTIN5boost7archive17archive_exceptionE` unmangles to:

```
typeinfo for boost::archive::archive_exception
```

What about the hex? Is it an address?

```
(gdb) x 0x0069ee50
0x69ee50 <_ZTVN10__cxxabiv121__vmi_class_type_infoE@@CXXABI_1.3+16>:    0x926be010
```

Yes, and the symbol unmangles to:

```
vtable for __cxxabiv1::__vmi_class_type_info
```



* * *

# tl;dr
* Life feels different now, I am on the other side. I have seen too much. I would very much like to never do this again.


[^8]: Thanks to [Dan](http://bytecove.co.uk/) for showing me [this symbol demangling tool](http://demangler.com/).

[^9]: A difficult task in itself as combinations of `assembly` `instruction` and `which way around` yielded many
    results containing questions pertaining to construction of flatpack furniture from a certain blue and yellow warehouse.
