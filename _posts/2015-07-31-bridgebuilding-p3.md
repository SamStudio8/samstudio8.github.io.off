---
layout: post
title: "The Tolls of Bridge Building: Part III, Sample (Un)Improvement"
---

[Previously, on Samposium]({% post_url page.previous.url %}): I finally had the **870 lanelets** required for the
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
%.coordinate_sort.bam: %.bam fixcoordsort$
  ${SAMTOOLS_BIN} sort -@ ${LSF_CPU} -n -T ${TMP_DIR}$(word 1,$+) -O bam -o $@ $<
[...]
```

`samtools sort` accepts an `-n` flag, to sort by query name, rather than the default co-ordinate position
(*i.e.* chromosome and base position on the reference sequence). Yet *somehow* the `Makefile` had been configured
to use query name sorting for both. I knew I was the last one to edit this part of the file, as I'd altered
the `-T` argument to prevent the temporary file clobbering discovered in the [last episode]({% post_url page.previous.url %}).

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

The indexing still failed. The lanelet bridged BAMs were still unsorted.

What? I double checked that I hadn't accidentally checked the last set of unsorted files.
Sadly I am somewhat competent and I'd attempted to index the right files, they really were unsorted.



### Vanishing Read Groups

* * *
#tl;dr
* Make sure your pipeline has the correct permissions to do what it needs with your directory, idiot.
* There's a special place in hell reserved for those who know `git` and don't use it[^3].
* [Brunel naively assumes inputs are sorted](https://github.com/wtsi-hgi/bridgebuilder/issues/7)


[^1]: Unfortunately somebody towing a caravan decided to have their car burst in to flame on the southbound M11, this and several other incidents turned a boring five hour journey in to a nine hour tiresome ordeal.

[^2]: Somewhat like a glorified `sudo` for humgen projects.

[^3]: Honestly, it doesn't matter how trivial the file is, if it's going to be messed around with frequently, or is a pinnacle piece of code for orchestrating a pipeline, put it under version control. Nobody is going to judge you for trivial use of version control.
