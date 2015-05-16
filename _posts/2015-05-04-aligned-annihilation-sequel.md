---
layout: post
title: "Aligned Annihilation II: Dumpster Diving"
excerpt: "I tried to extract a single integer from a core dump and instead fell in to an abyss and learned how to be a computer."
---

<p class="message"><b>Warning</b></br>This is going to get a bit technical and probably quite boring.</br>
If you are a victim of a boost::archive::archive_exception, skip to <a href="#fin">my conclusion</a>.<p>

If you haven't already looked at what led me to this awful situation, check out what happened when
I [annihilated all my alignment]({% post_url 2015-05-01-aligned-annihilation %}) data.

# The Hard Way
I've managed to get `rapsearch` to generate two core dumps for inspection with `gdb` but I am struggling to
extract any useful information beyond what I already know. The `abort()` and `raise()` can clearly be
seen in the backtrace and prior to that, construction of a standard library exception. There appears to be
a brief scene missing between `rapsearch` calling out to `boost::archive::basic_binary_iarchive` and
the error being constructed, leaving us with those `??` frames.

{% gist 8c88fdd06cd250872d67 backtrace.gdb %}

Desperately attempting to avoid having to edit and recompile `rapsearch`, I began nosing around the core
in a similar fashion as one would poke a stick around in a dirty pond. At first I naively tried to explore
frame 3, treating `0x31926bcbd6` as "the exception" before realising the address was for a function.
If we translate ("unmangle"[^8]) the symbol we can guess it is responsible
for handling assignment of an `exception_ptr`:

```
std::__exception_ptr::exception_ptr::operator=(std::__exception_ptr::exception_ptr const&)
```

I'd have to try harder.

I started looking at the registers for the same frame, as the function accepts an exception pointer as a parameter
it should be stored in a register here. It was when I started reading about [x86 Calling Conventions](http://en.wikipedia.org/wiki/X86_calling_conventions#System_V_AMD64_ABI) that I realised I was
probably in over my head. But a victim to the sunken cost fallacy now, I had to continue.
I disassembled the frame:

{% gist 8c88fdd06cd250872d67 disassemble_f3.gdb %}

It is [my understanding](http://stackoverflow.com/questions/2535989/what-are-the-calling-conventions-for-unix-linux-system-calls-on-x86-64)
that the `%rdi` register contains the first parameter given to a function and so here I'd expect to see
some form of pointer address. I'll dump the values held in the registers at this frame too.

{% gist 8c88fdd06cd250872d67 registers_f3.gdb %}

Unfortunately `0xb96` doesn't appear to point to anything. I wondered whether it would be
worth trying to manually work through the assembly instructions instead.
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
Beep boop. Ok, so that's a start. `_ZNSt15__exception_ptr13exception_ptrC1ERKS0_` unmangles[^8] to:

```
std::__exception_ptr::exception_ptr::exception_ptr(std::__exception_ptr::exception_ptr const&)
```

A constructor! Expecting a reference to pointer as its first parameter. So what's left in `%rdi`
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
Just `typeinfo`, not an actual `archive_exception` instance. Weird. What about the hex?
Is it an address? Where does it go?

```
(gdb) x 0x0069ee50
0x69ee50 <_ZTVN10__cxxabiv121__vmi_class_type_infoE@@CXXABI_1.3+16>:    0x926be010
```

It is, and the symbol unmangles to:

```
vtable for __cxxabiv1::__vmi_class_type_info
```

Hm. Too far. I'm not interested in the contents of a `vtable`. The `archive_exception` is a virtual class
so this will be where it's function pointers are populated at runtime. We want the specific instance of the class
that is raised, we want that error code.

Incidentally, eagle-eyed viewers will note the result of my computation was already stored in `%rbp`.
Oddly I thought this is where the frame's base pointer should be and find it unclear why it instead
points to a `typeinfo` object. But as disclaimed, I have no idea what I'm doing here.

Lunch slipped by me as I tried endless combinations of address lookups, quickly getting lost in the
37GB core file. Let's look back up the stack.

```
#6  0x00000000004452e2 in boost::archive::basic_binary_iarchive<boost::archive::binary_iarchive>::init() ()
```

Frame 6 holds the actual call that leaves `rapsearch` in an error state. Disassembly clocks in at about
250 lines of instructions so I'll cut out some potentially interesting lines instead:

{% gist 8c88fdd06cd250872d67 disassemble_f6.abridged.gdb %}

`_ZN5boost7archive17archive_exceptionC1ENS1_14exception_codeEPKcS4_` takes my eye and rightfully so,
it demangles to:

```
boost::archive::archive_exception::archive_exception(boost::archive::archive_exception::exception_code, char const*, char const*)
```

Bingo. We've got where the `archive_exception` is constructed! The first parameter is the `exception_code`.
It appears the contents of `%r13` are moved to `%rdi` just prior to the call.

Yet inspecting the various registers, I still can't get out the exception code.
Dan then asked me if I was sure the symbols are definitely correct. Recalling the recent recompilation debacle,
which was complete with its own `boost` library mis-version misadventure and the `rapsearch` repository
containing a two year old version of `boost`, no, I couldn't be sure at all.

I recompiled `rapsearch` and tried to run `gdb` on the same core. The symbols were different.
Shortly after, the work day drew to a close[^d] and sick of feeling like I'd been manually
searching through a cow's core dump, I gave up.

<blockquote>it's ok at least we had fun right
<footer>— <a href="http://bytecove.co.uk/">Daniel Evans</a></footer>
</blockquote>

# The Easy Way
There must be an easier way, I thought, just before bed. The [prelude to this epic]({% post_url 2015-05-01-aligned-annihilation %}) introduced the error at hand:

```bash
terminate called after throwing an instance of 'boost::archive::archive_exception'
  what():  invalid signature
/cm/local/apps/sge/current/spool/node012/job_scripts/1442997: line 40: 31100 Aborted
  (core dumped) rapsearch -q $QUERY -d /ibers/ernie/groups/rumenISPG/Databases/2015_04-trembl-ec3/rapsearch/2015-04__uniprot__ec_3__tax_2-Bacteria__reviewed_no.rap -u 1 -z 8 -e 0.00001 > $OUTFILE
```

That post also linked to a [boost serialization archive exceptions manual entry](http://www.boost.org/doc/libs/1_37_0/libs/serialization/doc/exceptions.html) that Dan kindly located.
Listed within are various types of errors that can be raised, mapped by an enum:

```C
    [...]
    typedef enum {
        [...]
        invalid_signature,      // first line of archive does not contain
                                // expected string
        [...]
    } exception_code;
    exception_code code;
    [...]
```

Wait. `what()`? That looks familiar?

```bash
terminate called after throwing an instance of 'boost::archive::archive_exception'
  what():  invalid signature
```

Oh dear[^10]. I'd disregarded the `what()` error as it sounded confusing and mysterious and not related
to file handling. Yet it was trying to tell us the answer all along:

<blockquote>
<b>invalid_signature</b></br>
Archives are initiated with a known string. If this string is not found when the archive is opened, It is presumed that this file is not a valid archive and this exception is thrown.
</blockquote>

<a name="fin"></a>
# So, what's the verdict?
I had a hunch that this might have been caused by corruption of the temporary files that `rapsearch` creates
in the working directory. For each job, `rapsearch` creates temp files in the format `<outname>.tmp<N>`, where
`outname` is the basename of the output file and `N` is the index of the temporary file, starting at 0.

It's not uncommon to execute multiple jobs that share an output directory. Here, I was trying
to keep my data organised by storing the alignment hits for bacterial, archaeal and fungal associated
hydrolases on my limpet contigs in the same place.

However, `rapsearch` creates a `.m8` storing hits along with a somewhat esoteric `.aln` alignment file.
But one cannot prevent the latter file from being generated (*Update*: Or so I thought at the time, [see below](#suppress_aln)). I thought I'd try and be clever and find
a way around having to just delete any `.aln` files once the job had completed and found `rapsearch`
accepts a `-u` option:

```
-u int : stream one result through stdout [1: m8 result, 2: aln result, default: don't stream any result through stdout]
```

Great. I'll specify `-u 1` for `.m8` output only and redirect `stdout` to my `$OUTFILE`. Job done.
Except not. When using this stream option, `rapsearch` isn't writing to any files and so has no
value to prepend to the `.tmp<N>` temp suffix. So what happens?

Every job ends up sharing a `.tmp0` file. Which as you can imagine goes down pretty well.
As a job progresses, `rapsearch` heads off to disk to
[write to its trusty temp file](https://github.com/zhaoyanswill/RAPSearch2/blob/95c866e9b818b7b4b9648ef4e0810a33300c3432/Src/mergeUnit.cpp#L43)
only to discover the archive header has been tampered with, which is upsetting enough to throw an error.
Cue the **invalid_signature** on stage error. Mystery solved.

Blackout. Drop curtain.

## Update: A few hours
All six jobs, running in harmony:

```
-rw-r--r--  1 msn users    0 May  5 03:26 contigs.2015-04__uniprot__ec_3__tax_2157-Archaea__reviewed_no.rap.rap6.wip.tmp0
-rw-r--r--  1 msn users    0 May  5 03:31 contigs.2015-04__uniprot__ec_3__tax_2157-Archaea__reviewed_yes.rap.rap6.wip.tmp0
-rw-r--r--  1 msn users    0 May  5 03:27 contigs.2015-04__uniprot__ec_3__tax_2-Bacteria__reviewed_no.rap.rap6.wip.tmp0
-rw-r--r--  1 msn users    0 May  5 03:27 contigs.2015-04__uniprot__ec_3__tax_2-Bacteria__reviewed_yes.rap.rap6.wip.tmp0
-rw-r--r--  1 msn users    0 May  5 03:24 contigs.2015-04__uniprot__ec_3__tax_4751-Fungi__reviewed_no.rap.rap6.wip.tmp0
-rw-r--r--  1 msn users    0 May  5 03:32 contigs.2015-04__uniprot__ec_3__tax_4751-Fungi__reviewed_yes.rap.rap6.wip.tmp0
```

Note those all important unique temporary file names!

## Update: Bedtime
In case this ever affects anybody else, I've notified the developers by opening
[an issue](https://github.com/zhaoyanswill/RAPSearch2/issues/17) on the `rapsearch` Github.
Hooray for open sorcery!

<a name="suppress_aln"></a>
## Update: A few days and another manpage later
Turns out, one can suppress the `.aln` file after all. As demonstrated in the `rapsearch`
[usage examples](https://github.com/zhaoyanswill/RAPSearch2/blob/95c866e9b818b7b4b9648ef4e0810a33300c3432/readme#L80),
the `-b` option (help entry below) can be set to 0.

```
-b int    : number of database sequence to show alignments [default: 100]. If it's -1, all results will be shown.
```

This seems somewhat counter-intuitive to me and is simply not something I had thought
of trying. If anything I'd have expected `-b 0` to just create an empty `.aln` file!
Nevertheless, this is a blog about metagenomics, not user interface design, so let's
[get on with some science]({{ page.next.url }}).

* * *

# tl;dr
* <strike>Life feels different now, I took it too far. I learned more about registers, calling conventions and assembly and have seen too much. I would very much like to never do this again.</strike>[^11]
* `rapsearch` probably corrupts temporary files causing job failure for searches writing to the same output directory if you don't use the `-o` option.
* I am not a very good computer.
* Read errors. **Believe** errors.


[^8]: Thanks to [Dan](http://bytecove.co.uk/) for showing me [this symbol demangling tool](http://demangler.com/), as well as for putting up with hours of remote interrogation. You may increment your drinks counter.

[^9]: Made more difficult as combinations of `assembly` `instruction` and `which way around` yielded many
    search results for construction of IKEA flatpack furniture.

[^d]: Two hours ago.

[^10]: An understatement.

[^11]: I can't justify complaint about staring in to the abyss of disassembly because my lack of attention brought on my fate.
