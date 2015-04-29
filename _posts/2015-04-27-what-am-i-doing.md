---
layout: post
title: "What am I doing?"
---

A week ago I had a progress meeting with
[Amanda](http://www.aber.ac.uk/en/cs/staff-list/staff_profiles/?staff_id=afc) and
[Wayne](http://www.aber.ac.uk/en/cs/staff-list/staff_profiles/?login=waa2); who make up the
supervisory team for the computational face of my project.

As Wayne had been away from meetings for a few weeks, I began with a roundup of everything
that has been going [disasterously wrong](http://samstudio8.github.io/2015/04/17/exits/)[^1].
Progress on a functional analysis of the [limpet data](http://samstudio8.github.io/2015/04/21/the-story-so-far-p1/)
has been repeatedly hindered by a lack of resources on our cluster which is simply strugging with the
sheer size of the jobs I'm asking of it.

# The Cluster Conundrum
I've encountered two main issues with job size here:

* Jobs that are large because the inputs are large but few (*e.g.* assembling raw reads contained in a pair of 42GB files), or
* Jobs that are large because although the inputs are small (< 100MB), there are thousands of them (*e.g.* `BLAST`'ing large numbers of contigs against a sharded database[^2])

## Small-Big Jobs
The former is somewhat unavoidable, if `velvet` wants to consume 450GB of RAM for an assembly and we 
want an assembly specifically from `velvet` then it's a case of having to wait patiently for one of the larger
nodes to become free enough to schedule the job. Indeed we could look for other assemblers[^3] and evaluate
their bold claims regarding reduced resource usage over competitors but often when we've found a tool that
*just works*, we like to keep things that way -- especially if we want to be able to compare results of other
assemblies that must be manufactured in the same way.

Cluster jobs require resources to be requested up-front and guesstimating (even generously) can often
lead to a job being terminated for exceeding its allowance; wasting queue time (days) as
well as execution time (days or weeks) and leaving you with nothing to show[^4].
The problem is in asking for too much, you queue for a node longer, but when finally scheduled
you effectively block others from 
using resources for a significant time period and [I'll make you feel bad for it](http://samstudio8.github.io/2015/04/26/memblame/).

Really the only way to get around these constraints is to try and minimize the dataset you have in the first place.
For example for assemblies you can employ one or more of:

**Normalization**
: Count appearances of substrings of length *k* (**k-mers**) present in the raw reads, then discard corresponding reads in a fashion that retains the distribution of k-mers. Discarding data is clearly lossy, but the idea is that as the distribution of k-mers is represented in the same way but with fewer reads.

**Partitioning**
: Attempt to construct a graph of all k-mers present in the raw reads, then partition it in to a series of subgraphs based on connectivity. Corresponding reads from each partition can then be assembled seperately and potentially merged afterwards. Personally I've found this method a bit hit and miss so far but would like to have time to investigate more.

**Subsampling**
: Select a more manageable proportion of reads from your dataset at random and construct an assembly. Not only very lossy, this in itself raises some interesting sampling bias issues (to go with your original environment sampling and PCR biases).

**Iterative Subsampling**
: Assemble a subsample from your data set and then align the contigs back to the original raw reads. Re-subsample from all remaining unaligned reads and create a second assembly, repeat the process until you have *N* different assemblies and are satisfied with the overall alignment (*i.e.* the set of remaining unaligned reads is sufficiently small). Tom in our lab group has been pioneering this approach and might hopefully give a better explanation of this than I can.

## Big-Small Jobs
The latter category is a problem actually introduced by trying to optimise cluster scheduling in the first place.
For example, an assembly can produce thousands of **contigs** ([groups of reads believed by an assembler to belong together](http://samstudio8.github.io/2015/04/23/the-story-so-far-p3/)) and often we want to know if
any interesting known sequences can be found on these contigs. Databases of interesting known sequences
are often [(very) large](http://samstudio8.github.io/2015/04/24/trembling/) and so to avoid submitting an inefficient long-running memory-hogging small-big job to locate thousands of different needles in thousands of different haystacks (*i.e.* `BLAST`'ing many contigs against a large database),
we can instead attempt to minimize the size of the job by amortising the work over many significantly smaller jobs.

For the purpose of `BLAST`[^6], we can shard both the contigs and the database of interesting sequences in to smaller pieces. This reduces the search space (fewer interesting-sequence needles to find in fewer contig haystacks) and thus execution time and resource requirements. Now your monolith job is represented by hundreds
(or thousands) of smaller, less resource intensive jobs that finish more quickly. Hooray!

Until the number of jobs you have starts [causing trouble](http://samstudio8.github.io/2015/04/17/exits/).

Of course this in turn makes handling data for downstream analysis a little more complex, output files need
converting, sorting and merging before potentially having to be re-sharded once again to fit them through a
different tool.

# Conquering Computational Complexity
So how can we move forward? 

# Ignoring the Impossible

* * *

# tl;dr
* My data is too big and my computer is too small.
* There are *big-small* jobs and *small-big* jobs and both are problematic and unavoidable.
* The project isn't impossible.

[^1]: Which is pretty much anything that involves a computer.

[^2]: In an attempt to speed up BLAST queries against large databases we have taken to splitting
    the database into 'shards'; submitting a job for each set of contigs against a specific
    database shard, before `cat`'ing all the results together at the end.
    
[^3]: In fact, currently I'm trying to evaluate [`MegaHIT`](https://github.com/voutcn/megahit).

[^4]: This isn't always strictly true. For example, aligners can flush output hits to a file as
   they go along and with a bit of fiddling you can pick up where you left off and `cat` the
   outputs together[^5].

[^5]: I call this **re-tailing**.

[^6]: Other short-read sequencer aligners are available.
