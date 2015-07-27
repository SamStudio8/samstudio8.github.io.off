---
layout: post
title: "The Tolls of Bridge Building: Part II, Construction [WIP]"
---

[Last time]({{ page.previous.url}}) on Samposium, I gave a more detailed look at the project I'm working
on and an overview of what has been done so far. We have **870** **lanelets** to pre-process and improve
into **samples**. In this post, I explain how the project has turned into a dangerous construction site.

While trying to anticipate expected work to be done in a [previous post]({% post_url 2015-06-24-sanger-sequel %}), I wrote:
> I also recall having some major trouble with needing to re-align the failed samples to a different reference: these
samples having failed, were not subjected to all the processing of their QC approved counterparts, which we'll need
to apply ourselves manually, presumably painstakingly.

'Painstakingly' could probably be described as an understatement now. For the past few weeks I've spent
most of my time overseeing the construction of bridges.

## Bridgebuilding
### References
Pretty much one of the first stops for genomic data spewed out of the instruments here is an alignment
to a known human reference. It's this step that allows most of our downstream analysis to happen, we
can "pileup" sequenced reads where they belong along the human genome, compare them to each-other
or other samples and attempt to make inferences about where variation lies and what its effects may be.

Without a reference -- like the metagenomic work I do back in Aberystwyth -- one has to turn to more
complicated de novo assembly methods and look for how the reads themselves overlap and line up in
respect to eachother rather than a reference.

Releases of the "canonical" human reference genome are managed by the
[Genome Reference Consortium](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/), of which
the Wellcome Trust Sanger Institute is a member. The GRC released `GRCh37` in March 2009, seemingly
also known as `hg19`. It is to this reference that our data was aligned to. However, to combat some
perceived issues with `GRCh37`, the 1000 Genomes Project team constructed their own derivative reference:
`hs37d5`.

To improve accuracy and reduce mis-alignment, `hs37d5` includes additional subsequences termed "decoy
sequences". These decoys prevent reads whose true origin is not adequately represented in the reference
from mapping to an arbitrarily similar looking location somewhere else in the genome. Sequences
of potential bacterial or viral contaminants, such as Epstein-Barr virus which is used in immortalisation
of some types of cell line are also included to prevent contaminated reads incorrectly aligning to the
human genome.

The benefits are obvious, fewer false positives and more mapped sequence for the cost of a re-alignment?
Understandably, our study decided to remap all good lanelets from `GRCh37` to `hs37d5` and re-run the
sample improvement process.

### Recipe
Of course, to execute our quality control pipeline fairly, we too must remap all of our lanelets to `hs37d5`.
Luckily, Josh and his team had already conuured some software to solve this problem in the form of
[`bridgebuilder`](https://github.com/wtsi-hgi/bridgebuilder) whose components and workflow are
modelled in the diagram below:

![]({{ site.url }}/public/posts/bridgebuilding-p2/bridge_builder_v1.png)

`bridgebuilder` has three main components summarised thus:

* **baker**  
  Responsible for the generation of the "reference bridge", mapping the old reference to
  the new reference. In our case, the actual `GRCh37` reference itself is bridged to `hs37d5`.
  The result is metaphorically a collection of bridges representing parts of `GRCh37` and their
  new destination somewhere in `hs37d5`.
* **binnie**  
  Takes the original lanelet as aligned to the old reference and using the `baker`'s reference
  bridge works out whether reads can stay where they are, or whether they fall in genomic regions
  covered by a bridge. For each of the 870 lanelets, the result is the population of four bins:
    * **Unbridged reads**  
      Reads whose location in `hs37d5` is the same as in `GRCh37` and do not require re-mapping.
      Unbridged reads are sub-binned by whether or not their ordering might have been changed.
    * **Bridged reads**  
      Reads that do fall on a bridge and are mapped across from `GRCh37` to `hs37d5`.
    * **Unmapped reads that are now mapped**  
      Reads that did not have an alignment to `GRCh37` but now align to `hs37d5`.
* **brunel**  
  Takes the `binnie` bins as "blueprints" and merges all reads to generate the new lanelet.

### Rinse and Repeat
Whilst designed to function together and co-existing in the same package, `bridgebuilder`
is not currently a *push button and acquire bridges* process. Luckily for me, back in 2014
while I was frantically finishing the write-up of my undergraduate thesis, Josh had **quickly**
prepared a **prototype** `Makefile` with the purpose of orchestrating our construction project.

This was neat, as it accepted a file of filenames (a "fofn") and automatically took care of
each step and all of the many intermediate files inbetween. More importantly it provided a
wrapper for handling submission of jobs to the compute cluster (the "farm"). The latter was
especially great news for me, having spent [most of my first year battling]({% post_url 2015-04-27-what-am-i-doing %}) our [job management system]({% post_url 2015-02-17-sun-grid-engine %}) back in Aberystwyth.

In theory remapping all 870 lanelets should have been accomplished through one simple
submission of `Make` to the farm. As you may have guessed by my tone, this was not the case.

### Ruins
Shortly before calling it a day back in 2014, Josh gave the orchestrating Makefile a few
shots. The results were hit and miss, a majority of the required files failed to generate
and the "*superlog*" containing over a million lines aggregated `stdout` from each of the submitted jobs
was a monolithic mish-mash of unnavigable segmentation faults, samtools errors, debugging
and scheduling information. Repeated application of the `Makefile` superjob yielded a handful
more desired files and the errors appeared to tail off, yet we were still missing several hundred
of the final lanelets. Diagnosing where the errors were arising was particular problematic due to
two assumptions in the albeit rushed design of our `Makefile`:

* **Error Handling**  
  The `Makefile` did not expect to have to handle the potential for error between one step
  and the next. Many but by no means all of the errors pumped to the superlog pertained to
  the absence of intermediate files, rather than a problem with execution of the commands themselves.
* **Error Tracing**  
  There was no convenient mechanism to trace back entries from the superlog to the actual
  cluster job (and thus particular lanelet) to try and reproduce and investigate errors.
  The superlog was just a means of establishing a heartbeat on the superjob and to roughly
  guage progress.

This is where the project left off and to where we must now return.

### Revisiting Ruins

...A major difference between then-and-now is the introduction of the `samtools quickcheck`
subcommand, written by Josh after the expected behaviour of `samtools index` was changed.
We could run all our current SAM/BAM alignment files through `quickcheck` to rule out basic
errors...





* often did not propagate errors with pipes
* stalled frequently
* assumed file creation as success


We had 870
**lanelets** -- parts of whole samples -- 


* Jobs not having enough time or memory
* Jobs failing to propagate an exit code
* Jobs failing stochastically
* 33 jobs failing as they had already been remapped
* 387 jobs failing as they were sat on top of an invalid reference

###
* missing RG (my fault)
* samtools sort clobbering bug
* "the 33"
* "177" files... empty bridge BAM
* 10% remainder... but...
* "2014" files, sai indexes built in 2014...
* "the 6"


* the final checks: quickcheck, 177, 2014, -e, -d UNKN, gzip -tf, samtools -c
* the final final checks: random subsample... all but bb-bam same size...
