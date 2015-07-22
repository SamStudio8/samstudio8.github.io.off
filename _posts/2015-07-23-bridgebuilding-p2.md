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
To fairly execute our quality control pipeline, when creating new improved samples for our fail-carrying
samples, we must follow all pre-processing steps from the prime study proper. Thus, we too must remap all of
our lanelets to `hs37d5`.

Luckily, Josh's team had in fact already solved this problem and created
[`bridgebuilder`](https://github.com/wtsi-hgi/bridgebuilder)

...it's components and workflow is modelled in the diagram below:
..![]({{ site.url }}/public/posts/bridgebuilding-p2/bridge_builder_v1.png)

* **baker**  
  ...generates the rules for bridging alignments given the origin and destination references...
* **binnie**  
  ...divides reads in to appropriate bins
* **brunel**  
  ...takes the "blueprints" and generates the final output...


### Rinse and Repeat
Back in mid-2014 while I was frantically finishing the write-up of my thesis, unknown to me
Josh prepared a `Makefile` with the purpose of orchestrating this process.

* did not attempt to specially detect errors
* often did not propagate errors with pipes
* did not print stderr
* log file could not determine which stdout came from which job
* stalled frequently
* assumed file creation as success

I'd hoped this had been mostly sorted before my return, but 

Unknown to me, back 

We had 870
**lanelets** -- parts of whole samples -- 


* Jobs not having enough time or memory
* Jobs failing to propagate an exit code
* Jobs failing stochastically
* 33 jobs failing as they had already been remapped
* 387 jobs failing as they were sat on top of an invalid reference

###
* missing RG
* "the 33"
* "177" files... empty bridge BAM
* 10% remainder... but...
* "2014" files, sai indexes built in 2014...
* "the 6"


* the final checks: quickcheck, 177, 2014, -e, -d UNKN, gzip -tf
* the final final checks: random subsample... all but bb-bam same size...
