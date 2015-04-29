---
layout: post
title: "What am I doing?"
---

A week ago I had a progress meeting with
[Amanda](http://www.aber.ac.uk/en/cs/staff-list/staff_profiles/?staff_id=afc) and
[Wayne](http://www.aber.ac.uk/en/cs/staff-list/staff_profiles/?login=waa2); who make up the
supervisory team for the computational face of my project, I talked about how computers are
terrible and where the project is heading.

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
[converting, sorting and merging](http://samstudio8.github.io/2015/04/29/pipelines/) before potentially having to be re-sharded once again to fit them through a different tool.

# Conquering Complications
So how can we move forward? We could just do what is [fashionable at the moment](https://xkcd.com/927/) and write 
a fantastic new [assembler|aligner|<software>] that is better and faster[^7] than the competition, uses next-to-no
memory and can even be run on a Raspberry Pi, but this is more than a PhD in itself[^8], so sadly, I guess
we have to make do and stick with what we have and attempt to use it more efficiently.

Digressing, I feel a major problem in bioinformatics software right now is a failure to adequately communicate
the uses and effects of parameters; how can end-users of your software fine tune[^9] controls and options
without it feeling like piloting a Soyuz? I think if the learning curve is too great, with understanding
hampered further by a lack of tutorials or extensive documentation with examples, users end up
driven to roll their own solution. Often in these cases the end result is maintained by a single developer
or group, missing out on the benefits of input from an open-source community at large.

Small-Big jobs can currently be tackled with novel methods like Tom's iterative subsampling
as described above, or of course, by adding additional resources (but that costs money).

Some of the risk recently identified with the execution of Big-Small jobs can be reduced by being a little
more organised. I'm in the process of writing some software to ease interaction with Sun Grid Engine that
now places logs generated during job execution outside of the working directory -- reducing some of the I/O load
when repeatedly requesting the contents of output directories.

Keeping abreast of the work of others who dared to tread and write their own new assembler, aligner or whatever
is important too. Currently we're testing out [`rapsearch`](http://omics.informatics.indiana.edu/mg/RAPSearch2/)
as an alternative to `BLAST` simply due to its execution speed (yet another post in itself).
`BLAST` is pretty old and "better" alternatives are known to exist, but it's still oft-cited and an expected part
of analysis in journal papers, so switching out parts of our pipeline for performance is not ideal. At the same
time, I actually want to get some work done and right now using `BLAST` on the dataset I have, with the resources
I have is proving too problematic.

At the very least, we can now use `rapsearch` to quickly look for hits to be analysed further with `BLAST`
if we fear that the community may be put off by our use of "non-standard" software.

# Ignoring the Impossible
After trading some graph theory with Wayne in return for some biological terminology, we turned our attention
to a broad view of where the project as whole is heading. We discussed how it is difficult to assemble entire
genomes from metagenomic datasets due to environmental bias, PCR bias and clearly, computational troubles.

I'd described my project at a talk previously;

<blockquote>...it's like trying to simultaneously assemble thousands of jigsaws but some of the jigsaws are heavily duplicated and some of the jigsaws hardly appear at all, a lot of the pieces are missing and quite a few pieces that really should fit together are broken. Also the jigsaws are pictures of sky.</blockquote>

Lately I've started to wonder how this is even possible, how can we state with confidence that we've
assembled a whole environment? How do we know the initial sample contained all the species? How can
we determine what is sequencing error and what is real and rare? How on Earth are we supposed to identify
all affinities in variation for all species across millions of reads that are shorter than my average
[Tweet](https://twitter.com/samstudio8/) that barely overlap?

We can't[^10].

But that's ok. That isn't the project. These sorts of aims are too broad, though that won't prevent
me from trying.
Currently I'm hunting for **hydrolases** (enzymes used to
break apart chemical bonds in presence of water), so we can turn the problem on its head a little.
Instead of creating an assembly and assigning taxonomic and functional annotations to every single one
of the resulting contigs then filtering the results by those that resemble hydrolasic behaviour
-- treating each contig as equally interesting -- we can just look at contigs that contain coding
regions for the creation of hydrolases directly! We can use a short-read aligner such as `rapsearch`
or `BLAST` to search for needles from a hydrolase-specific database of our own construction, instead
of a larger, more general bacterial database.

We can then query the assembly for the original raw reads that built the contig on which
strong hits for hydrolases appear. We can then take a closer look at these reads alone, filtering
out whole swathes of the assembly (and thus millions of reads) that are "not interesting" in terms
of our search.

We want to identify and extract interesting enzymes and the sequences that derive them, discovering a
novel species in the process is a nice bonus but the protein sequence is the key.

* * *

# tl;dr
* My data is too big and my computer is too small.
* There are *big-small* jobs and *small-big* jobs and both are problematic and unavoidable.
* There just isn't time to look at everything that is interesting.
* We need to know the tools we are using inside out and have a very good reason to make our own.
* We don't have to care about data that we aren't interested in.
* The project **probably** isn't impossible.

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

[^7]: Bonus points for ensuring it is also harder and stronger.

[^8]: I learned from my undergraduate dissertation that no matter how hard you try, the time to investigate every interesting side-street simply does not exist and it's important to try and stay on some form of track.

[^9]: I had a brief discussion about the difficulty of automated parameter selection on Twitter after a virtual conference and this is something I'd like to write more about at length in future:
    <blockquote class="twitter-tweet" lang="en"><p>.<a href="https://twitter.com/bioinformatics">@bioinformatics</a> nucleotid.es for biologists &amp; not to document &quot;fine tuning&quot;. But most assemblers need this tuning? <a href="https://twitter.com/hashtag/BaltiAndBioinformatics?src=hash">#BaltiAndBioinformatics</a></p>&mdash; Sam Nicholls (@samstudio8) <a href="https://twitter.com/samstudio8/status/557972285014155267">January 21, 2015</a></blockquote>
    <script async src="//platform.twitter.com/widgets.js" charset="utf-8"></script>

[^10]: Probably.
