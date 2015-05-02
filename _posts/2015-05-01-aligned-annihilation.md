---
layout: post
title: "Aligned Annihilation"
---

This afternoon in a coffee fueled fugue, I nuked every directory containing output for any attempt
to align the limpet contigs to any form of database so far[^1]. I was too far down the rabbit hole of which outputs
belonged to what job, which sets of outputs were supposed to be appended to the ends of other truncated
outputs, whether we were [looking at hydrolases or taxonomy]({% post_url 2015-04-27-what-am-i-doing %})
and were the hits to TrEMBL or SwissProt (or both, or neither)? Not that the database matters
given my copies of both are now [heavily redundant]({% post_url 2015-04-24-trembling %}) anyway.

As you might imagine, this was a briefly cathartic experience interrupted by an intense realisation
that I am definitely behind in producing interesting data. But it is refreshing to have
a grasp on what we have so far, even if it is now very little. It was also easy to work out where to
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

Of the two successful jobs, the only link appears to be their small database size (see table). Yet the smallest database
still met the same fateful error regardless, thus it appears the behaviour is intermittent.
I can only guess the trouble is related to the particular compute node rather than a bug in `rapsearch` or
an error in my job submission scripts (which have had at least six PhD candidate grade eyes on them now)
and so I've escalated my confusion and ruined the bank holiday weekend of our trusty sysadmin.

* * *

# tl;dr
* I deleted everything that hasn't quite worked and it felt good for a short time.
* The cluster is still not working.

[^1]: Apart from my first and only successful `BLAST` job that completed over Christmas that was against
    a heavily outdated SwissProt. I'm not even sure if this is of use any more but it is all I have to
    show for myself.

[^2]: Is that right?

[^3]: Especially useful for quickly re-submitting the jobs I sent with a missing database file extension.

[^4]: At least the core wasn't dumped this time. This also shows that my no error tolerance
    [configuration]({% post_url 2015-04-17-exits %}) works correctly, we caught the failure!
    
[^5]: I can't help but name drop my tool here. It just makes my life so much easier. Instead of messing
    around editing job submission scripts to change queues and parameters, I can just resubmit by answering
    a few prompts from the comfort of my terminal:

    ```bash
    sunblock resub 8-11
    ```
    
[^6]: Tom and I can't believe it either. It's the most free space I've seen since I started here:

    ```
    Filesystem                                           Size  Used Avail Use% Mounted on
    storage01.ib.cluster:/ibers/ernie/scratch             11T  6.5T  3.3T  67% /ibers/ernie/scratch
    ```
