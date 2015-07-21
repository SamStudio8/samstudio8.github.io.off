---
layout: post
title: "Securing your Six"
excerpt: "Sam finally realises his `apache2` stack is not working over IPv6."
---

As a financially constrained student, like many others, I use `apache`'s
[support for Server Name Indication (SNI)](https://wiki.apache.org/httpd/NameBasedSSLVHostsWithSNI)
to serve multiple SSL domains from one IP. I'm somewhat competent and the setup seems to work for
all of my domains. Yet, some time ago I tried to access one of my `VirtualHosts` from work over SSL
and was greeted by a fairly standard "invalid certificate" error. A certificate was produced but not for
the correct domain.

I had caused this by accident once before, where during a rushed deployment of new SSL keys following
Heartbleed, I was literally serving the wrong `SSLCertificateFile` to clients for that particular `VirtualHost`.
But after triple checking the configuration stanza, everything seemed to be correct in this instance.
What's more is the site was definitely receiving traffic and I could access it outside of work without error.

I dismissed the problem as a quirk of Sanger's network which has been known to do funny
things with web cache in the past, until a bug report from Germany rolled in. The same
website was not accessible from their home ISP on the continent.

"It works for me", I thought, and clearly for the majority of other users too. I could access
my other SSL protected domains from both work and so too could our bug reporting German counterpart.
I put it down to some weird quirk of Germany.

A little while later, I found that this particular website was still unaccessible from work.
Forcing a security exemption, the content downloaded is for the domain the certificate is for[^2].
Now far beyond any reasonable cache time, I figured something must really be wrong.

I scoured access and error logs, trying to find something obvious. I focused on the peculiar nature
of how other SSL protected domains worked fine and yet this one did not. I altered the `LogFormat` to
dump more information and finally noticed a discernable difference.

*The certificate error only occured when the client had an IPv6 address.*

Bollocks. I'd dun goofed the IPv6 configuration, pretty badly. Whilst the server itself is
"IPv6 ready": it can be pinged and the world's DNS servers know how to reach it over
the protocol, I'd never told `apache` it is expected to be able to serve content over SSL over IPv6.

After all the investigatory effort, the fix just consisted of a minor update the `apache` ports configuration
to add a new `NameVirtualHost` directive for the server's IPv6 address on both ports 80 and 443:

```
[...]
NameVirtualHost [<IPv6 Address>]:80
[...]

<IfModule mod_ssl.c>
    [...]
    NameVirtualHost [<IPv6 Address>]:443
</IfModule>
```

...and also to add the IPv6 address[^1] alongside the IPv4 to each of the `VirtualHost`:

```
<VirtualHost <IPv4 Address>:443 [<IPv6 Address>]:443>
    # Some configuration...
</VirtualHost>
```
**It works!**

* * *

#tl;dr
* I forgot to configure `apache` to serve content over IPv6 for SSL traffic, things went wrong.
* I configured `apache` with a `NameVirtualHost` for the server's IPv6 address and things are no longer wrong.

[^1]: Those square brackets aren't to be interpreted as "optional", they are how `apache` expects an IPv6 address to be formatted.

[^2]: The typical behaviour of `apache` not knowing which `VirtualHost` is supposed to be responding to a request is loading the first one.
