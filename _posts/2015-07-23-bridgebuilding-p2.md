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

## Bridgebuilding
Pretty much one of the first stops for genomic data spewed out of the instruments here is an alignment
to a known human reference. It's this step that allows most of our downstream analysis to happen, we
can "pileup" sequenced reads where they belong along the human genome, compare them to each-other
or other samples and attempt to make inferences about where variation lies and what its effects may be.

Without a reference -- like the metagenomic work I do back in Aberystwyth -- one has to turn to more
complicated de novo assembly methods and look for how the reads themselves overlap and line up in
respect to eachother rather than a reference.

Releases of the canonical human reference genome are managed by the
[Genome Reference Consortium](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/), of which
the Wellcome Trust Sanger Institute is a member.

There is however, more than one human reference sequence
...during the course of the study from which our data is sourced from, the decision was made to switch
from `GRCh37` to `hs37d5`



I'd hoped this had been mostly sorted before my return, but 
Josh did not have the time to
keep pumping 
'Painstakingly' could probably be described as an understatement now. For the past few weeks I've spent
most of my time overseeing the construction of bridges.



Unknown to me, back in 2014 while I was frantically finishing the write-up of my thesis, Josh prepared
a `Makefile` with the purpose of carrying out one of the more difficult tasks of the project. 

We had 870
**lanelets** -- parts of whole samples -- 


* Jobs not having enough time or memory
* Jobs failing to propagate an exit code
* Jobs failing stochastically
* 33 jobs failing as they had already been remapped
* 387 jobs failing as they were sat on top of an invalid reference
