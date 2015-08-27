---
layout: post
title: "The Tolls of Bridge Building: Part III, Sample (Un)Improvement"
---

[Previously, on Samposium]({{ page.previous.url }}): I finally had the **870 lanelets** required for the
sample improvement process. But in this post, I explain how my deep-seated paranoia in the quality of my
data just wasn't enough to prevent what happened next.

I submitted my 870 bridged BAMs to `vr-pipe`, happy to essentially be rid of having to deal with the data
for a while. `vr-pipe` is a complex bioinformatics pipeline that in fact consists of a series of pipelines
that each perform some task or another with many steps required for each. The end result is "improved sample"
BAMs, though perhaps due to the nature of our inclusion of failed lanelets we should title them
"unimproved sample" BAMs.
Having just defeated the supposedly automated bridging process was pretty happy to not have to be doing
this stuff manually and could finally get on with something else for a while... Or so I thought.

### A Quiet Weekend
It was Friday and I was determined to leave on-time to see my family in Swansea at a reasonable hour[^1].
After knocking up a quick Python script to automatically generate the metadata file `vr-pipe` requires,
describing the lanelets to be submitted and what sample they should be included in the "improvement" of,
I turfed the job over to Micheal who has fabled `mercury` access[^2] and could set-up `vr-pipe` for me.

Being the sad soul that I am, I occasionally checked in on my job over the weekend via the `vr-pipe` web interface,
only to be confused by the apparent lack of progress. The pipeline appeared to just be juggling jobs around
various states of pending. But without `mercury` access to inspect more, I was left to merely enjoy my weekend.

Between trouble on the server farm and the fact that my job is not particularly high priority,
the team suggested I be more patient and so I gave it a few more days before pestering Josh to take a look.

### Pausing for Permissions
As suspected, something had indeed gone wrong. Instead of telling anybody, `vr-pipe` sat on a small mountain
of errors, hoping the problem would just go away. I've come to understand this is expected behaviour.
Delving deep in to the logs revealed the simple problem: `vr-pipe` did not have sufficient write permissions
to the directory I had provided the lanelet files in, because I didn't provide group-write access to it.

One ```chmod 775 .``` later and the pipeline burst in to life, albeit very briefly, before painting the `vr-pipe`
web interface bright red. Evidently, the problem was more serious.

### Sorting Names and Numbers
The first proper step for `vr-pipe` is creating an index of each input file with `samtools index`.
Probably the most important thing to note for an index file, is that to create one, your file must be
correctly sorted. Micheal checked the logs for me and found that the indexing job had failed on all but
33 (shocker) of the input files, complaining that they were all unsorted.

But how could this be? `brunel` requires sorted input to work, our orchestrating `Makefile` takes
care of sorting files as and when needed with `samtools sort`. There must be some mistake!

I manually invoked `samtools index` on a handful of my previous bridged BAMs and indeed, they are unsorted.
I traced back through the various intermediate files to see when the sort was damaged before finally
referring to the `Makefile`. My heart sank:

```Make
[...]
# sort by queryname
%.queryname_sort.bam: LSF_MEM=10000
%.queryname_sort.bam: LSF_CPU=4
%.queryname_sort.bam: %.bam
  ${SAMTOOLS_BIN} sort -@ ${LSF_CPU} -n -T ${TMP_DIR}$(word 1,$+) -O bam -o $@ $<

# sort by coordinate position
%.coordinate_sort.bam: LSF_MEM=10000
%.coordinate_sort.bam: LSF_CPU=4
%.coordinate_sort.bam: %.bam
  ${SAMTOOLS_BIN} sort -@ ${LSF_CPU} -n -T ${TMP_DIR}$(word 1,$+) -O bam -o $@ $<
[...]
```

`samtools sort` accepts an `-n` flag, to sort by query name, rather than the default co-ordinate position
(*i.e.* chromosome and base position on the reference sequence). Yet *somehow* the `Makefile` had been configured
to use query name sorting for both. I knew I was the last one to edit this part of the file, as I'd altered
the `-T` argument to prevent the temporary file clobbering discovered in the [last episode]({{ page.previous.url }}).

Had I been lazy and naughty and copied the line from one stanza to the next? I was sure I hadn't.
But the knowing-grin Josh shot at me when I showed him the file, had me determined to try and prove my innocence.
Although, it should be first noted that a spot in a special part of hell must be first reserved for the both of us,
as the `Makefile` was not under version control.

I'd squirrelled away many of the original files from 2014, including the intermediates and selecting any
co-ordinate sorted file yielded a set of query name sorted reads. My name was cleared! Of course, whilst
this was apparently Josh's mistake, it's not entirely fair to point the finger given I never noticed the bug
despite spending more time nosing around the `Makefile` than anyone else. As mentioned, I'd even obliviously edited right next to the extraneous `-n` in question.

But I am compelled to point the finger elsewhere: `brunel` **requires** sorted inputs to work correctly, else it
creates files that clearly can't be used in the sample (un)improvement pipeline! How was this allowed to happen?

Sadly, it's a matter of design. `brunel` never anticipated that somebody might attempt to provide incorrect
input and just gets on with its job regardless. Frustratingly, had this been a feature of `brunel`, we'd have
caught this problem a year ago on the first run of the pipeline.

The co-ordinate sorting step does at least immediately precede invocation of `brunel`, which is relatively
fast. So after correcting the Makefile, nuking the incorrectly sorted BAMs and restarting the pipeline, it
wasn't a long wait before I was ready to hassle somebody with `mercury` access to push the button.

### Untranslated Translations
Before blindly resubmitting everything, I figured I'd save some time and face by adding
`samtools index` to the rest of my checking procedures, to be sure that this first indexing
step would at least work on `vr-pipe`.

The indexing still failed. The final lanelet bridged BAMs were still unsorted.

Despite feeding our now correctly co-ordinate sorted inputs to `brunel`, we still get an
incorrectly sorted output -- clearly a faux pas with `brunel`'s handling of the inputs.
Taking just one of the 837 failed lanelets (all but the already mapped 33 lanelets failed
to index) under my wing, I ran `brunel` manually to try and diagnose the problem.

The error is a little different from the last run, whereas before `samtools index` complained
about co-ordinates appearing out of order, this time, chromosomes appear non-continuously.
Unfortunately, several `samtools` errors still do not give specific information and this is one
of them. I knew that *somewhere* a read, on *some* chromosome appears somewhere it shouldn't, but
with each of these bridged BAMs containing tens of millions of records, finding it manually could be
almost impossible.

I manually inspected the first 20 or so records in the bridged BAM, all the reads were on TID 1,
as expected. The co-ordinates were indeed sorted, starting with the first read at position 10,000.

10,000? Seems a little high? On a hunch I grepped the bridged BAM to try and find a record on
TID 1 starting at position 0:

```
grep -Pn -m1 "^[A-z0-9:#]*\t[0-9]*\t1\t1\t" <(samtools view 7500_7#8.GRCh37-hs37d5_bb.bam)
```

I got a hit. Read `HS29_07500:7:1206:5383:182635#8` appears on TID 1, position 1. Its line
number in the bridged BAM? 64,070,629. There's our non-continuous chromosome.

I took the read name and checked the three sorted inputs to `brunel`. It appears in the "unchanged"
BAM. You may recall these "unchanged" reads are those that `binnie` deemed as not requiring re-mapping
to the new reference and can effectively "stay put". The interesting part of the hit? In the unchanged BAM,
it doesn't appear on chromosome 1 at all but on "HSCHRUN_RANDOM_CTG19", presumably some form of
decoy sequence.

This appears to be a serious bug in `brunel`. This "HSCHRUN_RANDOM_CTG19" sequence (and presumably others)
seem to be leaving `brunel` as aligned to chromosome 1 in the new reference. Adding some weight to the theory,
the unchanged BAM is the only input that goes through translation too.

Let's revisit [`build_translation_file`](https://github.com/wtsi-hgi/bridgebuilder/blob/673f8aa57abe7acbd688accbeda163c9f6b4eb6a/brunel/src/main.c#L86), you recall that the i-th SQ line of the input BAM -- the "unchanged" BAM --
is mapped to the j-th entry of the user-provided translation table text file. The translation itself
is recorded with `trans[i] = j` where both `i` and `j` rely on two loops meeting particular exit conditions.

But note the while loop:

```C
int counter = file_entries;
[...]
while (!feof(trans_file) && !ferror(trans_file) && counter > 0) {
    getline(&linepointer, &read, trans_file);
    [...]
    trans[i] = j;
    counter--;
}
[...] 
```

This while loop, that houses the `i` and `j` loops, as well as the assignment of `trans[i] = j`,
may not run for as many SQ lines (`file_entries`) found in the unchanged BAM, as stored in counter,
if the length of the `trans_file` is shorter than the number of `file_entries`. Which in our case, it is
-- as there are only entries in the translation text file for chromosomes that actually need translating (Chr1 -> 1).

Thus not every element in `trans` is set with a value by the while loop. Though we've 
[already seen]({{ post.previous.url }}) where there is no translation to be made, this doesn't work either.

This is particularly troubling, as the
[check for whether a translation should be performed](https://github.com/wtsi-hgi/bridgebuilder/blob/673f8aa57abe7acbd688accbeda163c9f6b4eb6a/brunel/src/main.c#L279)
relies on the default value of `trans` being `NULL`.
As no default value is set and 0 is the most likely value to turn up in the memory
allocated by `malloc`, the default value for `trans[i]` for all `i`, is 0. Or in `brunel` terms, **if I
can't find a better translation, translate SQ line `i`, to the 0th line (first chromosome) in the user
translation file**.

Holy crap. [That's a nasty bug](https://github.com/wtsi-hgi/bridgebuilder/issues/8).

As described in the bug report, I toyed with quick fixes based on initialising `trans` with default values:

> * **Initialise `trans[i] = i`**  
> This will likely incorrectly label reads as aligning to the i-th SQ line in the new file. If the number of SQ lines in the input is greater than that of the output, this will also likely cause samtools index to error or segfault as it attempts to index a TID that is outside the range of SQ.
> * **Initialise `trans[i] = NULL`**  
> Prevents translation of TID but actually has the same effect as above
> * **Initialise `trans[i] = -1`**  
> The only quick-fix grade solution that works, causes any read on a TID that has no translation to be regarded as "unmapped". Its TID will be set to "*" and the read is placed at the end of the result file. The output file is however, valid and indexable.

In the end the question of where these reads should actually end up is a little confusing. Josh seems to think
that the new bridged BAM should contain all the old reads on their old decoy sequences if a better place in the
new reference could not be found for them. In my opinion, these reads effectively "don't map" to hs37d5 as they
still lie on old decoy sequences not found in the new reference, which is convinient as my `trans[i] = -1`
initialisation marks all such reads as unmapped whilst also remaining the most simple fix to the problem.

Either way, the argument as to what should be done with these reads is particularly moot for me and the
QC study, because we're only going to be calling SNPs and focussing on our
[**Goldilocks Region**]({% post_url 2015-07-22-bridgebuilding %}) on Chromosome 3, which is neither a decoy
region or Chromosome 1.

Having deployed my [own fix](https://github.com/SamStudio8/bridgebuilder/tree/fix6), I re-ran our pipeline
for what I hoped to be the final time, re-ran the various checks that I've picked up along the way (including
`samtools index`) and was finally left with **870 lanelets** ready for unimprovement.

### Vanishing Read Groups
Or so I thought.

Having convinced Micheal that I wouldn't demand he assume the role of `mercury` for me at short notice again,
he informed me that whilst all my lanelets had successfully passed the *first step* of the *first pipeline*
in the entire `vr-pipe` workflow -- indexing. All but (you guessed it), 33 lanelets, failed the following
step. `GATK` was upset that some reads were missing an `RG` tag.

Sigh. Sigh. Sigh. Table flip.

`GATK` at least gave me a read name to look up with grep and indeed, these reads *were* missing their RG tag.
I traced the reads backward through each intermediate file to see where these tags were lost. I found that
these reads had been binned by `binnie` as requiring re-mapping to the new reference with `bwa`. The resulting
file from `bwa` was missing an `@RG` line and thus each read had no `RG` tag.

Crumbs. I hit the web to see whether this was anticipated behaviour. The [short answer](http://gatkforums.broadinstitute.org/discussion/1903/question-about-rg-tags) was yes, but the longer
answer was "you should have used `-r` to provide `bwa` with an `RG` line to use". Though I throw my hands up
in the air a little here to say "I didn't write the `Makefile` and I've never used `bwa`".

Luckily, Martin has recently [drafted a pull request](https://github.com/samtools/samtools/pull/371) to
`samtools` for a new subcommand: `addreplacerg`. Its purpose? Adding and replacing RG lines and tags in
BAM files. More usefully, at least to us, <strike>its default</strike>[^4] it offers an operation "mode" to tag "orphan" records (reads
that have no RG line -- exactly the problem I am facing) with the first RG line found in the header.

Perfect. I'll just feed each final bridged BAM I have to `samtools addreplacerg` and maybe, just maybe,
we'll be able to pull the trigger on `vr-pipe` for the final time.

I hope.

<a name="queuequeue"></a>
### Queue Queue *Update: 1 day later*
For those still following at home, my run of `samtools addreplacerg` seemed to go without a hitch, which
in retrospect should have been suspicious. I manually inspected just a handful of the files to ensure both
the "orphaned" reads now had an `RG` tag (and more specifically that it was the correct and only `RG`
line in the file) and that the already tagged reads had not been interfered with. All seemed well.

After hitting the button on `vr-pipe` once more, it took a few hours for the web interface to catch up and
throw up glaring red progress bars. It seems the first step - the BAM indexing was now failing? I had
somehow managed to go a step backwards?

The bridged BAMs were truncated... Immediately I begun scouring the source code of `samtools addreplacerg`
before realising in my excitement I had skipped my usual quality control test suite. I consulted
`bhist -e`, a command to display recently submitted cluster jobs that had exited with error and was
bombarded with line after line of `addreplacerg` job metadata. Inquiring specifically, each and every
job in the array had violated its run time limit.

I anticipated `addreplacerg` would not require much processing, it just iterates over the input BAM,
slapping an `RG` sticker on any record missing one and throws it on the output pile. Of course, with
tens of millions of records per file even the quickest of operations can aggregate into considerable time.
Thus placing the `addreplacerg` jobs on Sanger's `short` queue was clearly an oversight of mine.

I resubmitted to the `normal` queue which permits a longer default run time limit and applied
quality control to **all** files. We had the green light to re-launch `vr-pipe`.

### Conquering Quotas *Update: 4 days later*
Jobs slowly crawled through the indexing and integrity checking stages and eventually began their way
through the more time-consuming and intensive GATK indel discovery and re-alignment steps, before
suddenly and uncontrollably failing in their hundreds. I watched a helpless `vr-pipe` attempt to
resuscitate each job three times before calling a time of death on my project and shutting down the pipeline
entirely.

Other than debugging output from `vr-pipe` scheduling and submitting the jobs, the error logs were empty.
`vr-pipe` had no idea what was happening, and neither did I. 

I escalated the situation to Martin, who could dig a little deeper with his `mercury` mask on. The problem
appeared somewhat more widespread than just my pipeline; in fact all pipelines writing to the same
storage cluster had come down with a case of sudden unexpected rapid job death. It was serious.

The situation: pipelines are orchestrated by `vr-pipe`, which is responsible for submitting jobs
to the LSF scheduler for execution. LSF requires a user to be named as the owner of a job and so `vr-pipe`
uses `mercury`. I am unsure whether this is just practical, or whether it is to ensure jobs get a fair share
of resources by all having the same owner though I suspect it could just be an inflexibility in `vr-pipe`.
The net result of jobs being run as `mercury` is that every output file is also owned by `mercury`.
The relevance of this is that every storage cluster has user quotas, an upper bound on the amount
of disk space files owned by that user may occupy before being denied write access.

Presumably you can see where this is going. In short, my job pushed `mercury` 28TB over-quota and so
disk write operations failed. Jobs, unable to write to disk aborted immediately but the nature of the
error was not propagated to `vr-pipe`.

Martin is kindly taking care of the necessary juggling to rebalance the books. Will keep you posted.


* * *
#tl;dr
* Make sure your pipeline has the correct permissions to do what it needs with your directory, idiot.
* There's a special place in hell reserved for those who know `git` and don't use it[^3].
* [Brunel naively assumes inputs are sorted](https://github.com/wtsi-hgi/bridgebuilder/issues/7)
* [[critical] Brunel unintentionally translates untranlated entries](https://github.com/wtsi-hgi/bridgebuilder/issues/8)
* You should use the `-r` argument of `bwa sampe` to ensure your resulting BAM has an RG line
* Bioinformatics is *never* **ever** simple.
* There are bugs **everywhere**, you probably just haven't found them yet.


[^1]: Unfortunately somebody towing a caravan decided to have their car burst in to flame on the southbound M11, this and several other incidents turned a boring five hour journey in to a nine hour tiresome ordeal.

[^2]: Somewhat like a glorified `sudo` for humgen projects.

[^3]: Honestly, it doesn't matter how trivial the file is, if it's going to be messed around with frequently, or is a pinnacle piece of code for orchestrating a pipeline, put it under version control. Nobody is going to judge you for trivial use of version control.

[^4]: [Turns out](https://github.com/mp15/samtools/pull/4), its default mode is `overwrite_all`.
