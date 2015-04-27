---
layout: post
title: "What am I doing? [WIP]"
---

<p class="message">I'll finish it soon I promise.</message>

A week ago I had a progress meeting with
[Amanda](http://www.aber.ac.uk/en/cs/staff-list/staff_profiles/?staff_id=afc) and
[Wayne](http://www.aber.ac.uk/en/cs/staff-list/staff_profiles/?login=waa2); who make up the
supervisory team for the computational face of my project.

# A Recap
As Wayne had been away from meetings for a few weeks, I began with a roundup of everything
that has been going [disasterously wrong](http://samstudio8.github.io/2015/04/17/exits/)[^1].
Progress on a functional analysis of the [limpet data](http://samstudio8.github.io/2015/04/21/the-story-so-far-p1/)
has been repeatedly hindered by a lack of resources on our cluster which is simply strugging with the
sheer size of the jobs I'm asking of it.

I've encountered two main issues with job size here:
* Jobs that are large because the inputs are large (*e.g.* assembling raw reads contained in a pair of 42GB files), or
* Jobs that are large because although the inputs are small (< 100MB), there are thousands of them (*e.g.* `BLAST`'ing large numbers of contigs against a sharded database[^2])

# The Cluster Conundrum
The former is somewhat unavoidable, if `velvet` wants to consume 450GB of RAM for an assembly and we 
want an assembly specifically from `velvet` then it's a case of having to wait patiently for one of the larger
nodes to become free enough to schedule the job. Indeed we could look for other assemblers[^3] and evaluate
their bold claims regarding reduced resource usage over competitors but often when we've found a tool that
*just works*, we like to keep things that way -- especially if we want to be able to compare results of other
assemblies that must be manufactured in the same way.

Cluster jobs require resources to be requested up-front and guesstimating (even generously) can often
lead to a job being terminated for exceeding its allowance; wasting queue time (days) as
well as execution time (days or weeks) and leaving you with nothing to show[^4].
The problem is in asking for too much, you queue for a node longer and effectively block others from 
using resources for a significant time period and [I'll make you feel bad for it](http://samstudio8.github.io/2015/04/26/memblame/).

* * *

[^1]: Which is pretty much anything that involves a computer.

[^2]: In an attempt to speed up BLAST queries against large databases we have taken to splitting
    the database into 'shards'; submitting a job for each set of contigs against a specific
    database shard, before `cat`'ing all the results together at the end.
    
[^3]: In fact, currently I'm trying to evaluate [`MegaHIT`](https://github.com/voutcn/megahit).

[^4]: This isn't always strictly true. For example, aligners can flush output hits to a file as
   they go along and with a bit of fiddling you can pick up where you left off and `cat` the
   outputs together[^5].

[^5]: I call this **re-tailing**.
