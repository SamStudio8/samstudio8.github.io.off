---
layout: post
title: "Aligned Annihilation"
---

This afternoon in a coffee fueled fugue, I nuked every directory containing output for any attempt
to align the limpet contigs to any form of database so far. Here's why, and what I did next.

I awoke in the morning knowing I was too far down the rabbit hole to truly know which outputs
belonged to what job, which sets of outputs were now supposed to be appended to the ends of other truncated
outputs, whether I was even supposed to be
[looking at hydrolases or taxonomy]({% post_url 2015-04-27-what-am-i-doing %})
and did I try hitting on TrEMBL or SwissProt (or both, or neither)? Not that the database matters
given my copies of both are now [heavily redundant]({% post_url 2015-04-24-trembling %}) anyway.

As you might imagine, the obliteration was a briefly cathartic experience interrupted too soon by
an intense realisation that I am definitely behind in producing interesting data. But it is refreshing to have
a grasp on what we have so far, even if it is now very little[^1]. It was also easy to work out where to
go from here: generate some damn data.

Recently I [sidelined `BLAST`]({% post_url 2015-04-27-what-am-i-doing %}) in favour of `rapsearch`
for performing searches over contigs for "interesting sequences" as the cluster simply appears
[unable to cope]({% post_url 2015-02-17-sun-grid-engine %}) with such scale: once both the limpet assembly contigs
had been sharded into 250 bite-size pieces and hit against say, 30 shards of TrEMBL v2015_03,
you're already looking at 7,500 jobs -- which as [seen previously]({% post_url 2015-04-17-exits %}) is
potentially problematic.

Having re-stocked my hydrolase-specific databases after the latest release of Uniprot, with bacterial
and archaeal[^2] hydrolases, Wayne suggested to grab those from fungi too -- as they are an expected
taxonomic group in the dataset and we'd like to not miss them out. The table below summarises the
databases for future reference:

| Taxonomy | Reviewed | #Records  | DB Size (MB) | Failed |
|----------|:--------:|-----------|--------------|:------:|
| Bacteria | N (Tr)   | 1,069,711 | 2,700~       | X      |
| Bacteria | Y (SP)   | 42,027    | 108          | X      |
| Archaea  | N (Tr)   | 20,456    | 59           |        |
| Archaea  | Y (SP)   | 2,132     | 20           | X      |
| Fungi    | N (Tr)   | 61,521    | 243          | X      |
| Fungi    | Y (SP)   | 5,546     | 37           |        |

Thanks to the speed of `rapsearch`, one does not need to perform any sharding on the contig queries
or the database, which is lovely. Additionally, my
[work in progress](https://github.com/samstudio8/sunblock): `Sunblock`,
takes some of the pain out of working with the cluster, prompting a few questions
before automatically generating the job script and submitting it to SGE[^3] for execution.

Unfortunately the smooth sailing stopped shortly after scheduling and four of the six jobs terminated given
a non-zero exit (as configured) from the `rapsearch` command and dumping the following [familiar library error]({% post_url 2015-04-17-exits %})[^4]:

```bash
terminate called after throwing an instance of 'boost::archive::archive_exception'
  what():  invalid signature
/cm/local/apps/sge/current/spool/node012/job_scripts/1442997: line 40: 31100 Aborted
  (core dumped) rapsearch -q $QUERY -d /ibers/ernie/groups/rumenISPG/Databases/2015_04-trembl-ec3/rapsearch/2015-04__uniprot__ec_3__tax_2-Bacteria__reviewed_no.rap -u 1 -z 8 -e 0.00001 > $OUTFILE
```

A quick search offers no concrete solution. The disk isn't full[^6] and the accounting file indicates that
both the RAM and time allocated to the jobs were sufficient. You may recall my encountering this `boost` error
recently when client nodes had been upgraded in an attempt to address the "recent" kernel NFS bug disaster.

`node012` is part of the `large.q` and so I thought perhaps this heavily utilised big-job queue was yet to receive
(or yet to reboot following) such an update. After resubmitting the jobs to other queues with `Sunblock`[^5]
and found the behaviour repeated on another node, my hunch was proven incorrect. Grr.

Of the two successful jobs, the only link appears to be their small database size (see table). Yet the smallest database still met the same fateful error regardless, thus it appears the behaviour is unpredictable.
I can only guess the trouble is related to a particular compute node rather than a bug in `rapsearch` or
an error in my job submission scripts (which have had at least six PhD candidate grade eyes on them now)
and so I've escalated my confusion and ruined the bank holiday weekend of our trusty sysadmin.

## Update: The following morning
I spoke to Tom about the `rapsearch` jobs he's running right now and all is apparently well
on `node009` and `node010`. So I figured I'd resubmit the previously failed fungi TrEMBL job
to `node010` and wait for a result: it completed sucessfully in just over 40 minutes.

Unhelpfully, the fungi SwissProt job (which admittedly was previously successful) was able to
finish on the suspicious `node003`. So the problem is unpredictable and transient? The plot thickens.

I'll try out the bacterial and failed archaea tasks on `node010`.

## Update: Mostly still morning
The archaeal task was a success on `node007`, taking just shy of 20 minutes. One of the bacterial
re-resubmissions has already failed with the same error as before. I decided to have a nose around the `rapsearch` [repository](https://github.com/zhaoyanswill/RAPSearch2/tree/95c866e9b818b7b4b9648ef4e0810a33300c3432)
to try and track down where calls to `boost` archive serialization functionality even reside.
Instead I'm sat here confused as to why there is a `boost/` directory included in the source that hasn't been
updated for two years, whilst the project's [Makefile](https://github.com/zhaoyanswill/RAPSearch2/blob/95c866e9b818b7b4b9648ef4e0810a33300c3432/Src/Makefile) still links against a local installation of `boost` with `-lboost_serialization` as a sane person would? Wat?
I'm not even sure what the consequences of this are, does this old directory get used during building? 
I'll have to ask our sysadmin how he managed to build it as a module but this does not sit well with me. 

My good friend and resident expert in all things C(++), [Dan](http://bytecove.co.uk/) located some manual pages
that I should have read regarding [boost serialization exceptions](http://www.boost.org/doc/libs/1_37_0/libs/serialization/doc/exceptions.html).
The `boost::archive_exception` object holds an `exception_code` that maps to an enum to tell you exactly
what went wrong, which is awesomely thought out and very useful. Unfortunately for me, the error isn't caught
by `rapsearch` which is why I have an upset `stderr`. My options are probably either to re-compile
`rapsearch` with my own attempt at error handling, or re-enabling core dumps and inspecting the steaming remains
with `gdb`. I'm opting for the latter because quite frankly, screw spending my weekend 
messing around [trying to diagnose and fix](https://github.com/samtools/samtools/pull/259) poorly documented
C-based bioinformatics tools again.

## Update: A few days later
I've managed to get `rapsearch` to generate two core dumps for inspection with `gdb` but I am struggling to
extract any useful information beyond what I already know. The backtrace:

{% gist 8c88fdd06cd250872d67 backtrace.gdb %}

* * *

# tl;dr
* I deleted everything that hasn't quite worked and it felt good for a short time.
* The cluster continues to fail in an unexpected and unpredictable fashion.
* If your library gives you nice exception objects, try to catch them and use them for good.

[^1]: I kept my first and only successful `BLAST` job which completed over Christmas, that was against
    a heavily outdated SwissProt. I doubt at this point it is even of any use but it is all I have to
    show for myself.

[^2]: Is that right?

[^3]: Especially useful for quickly re-submitting the jobs I sent with a missing database file extension.

[^4]: At least the core wasn't dumped this time. This also shows that my no error tolerance
    [configuration]({% post_url 2015-04-17-exits %}) works correctly, we caught the failure!
    
[^5]: I can't help but name drop my tool here. It just makes my life so much easier. Instead of messing
    around editing job submission scripts to change queues and parameters, I can just resubmit by answering
    a few prompts from the comfort of my terminal!

    ```bash
    sunblock resub 8-11
    ```
    
[^6]: Tom and I can't believe it either. It's the most free space I've seen since I started here:

    ```
    Filesystem                                           Size  Used Avail Use% Mounted on
    storage01.ib.cluster:/ibers/ernie/scratch             11T  6.5T  3.3T  67% /ibers/ernie/scratch
    ```
