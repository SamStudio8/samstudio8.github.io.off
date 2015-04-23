---
layout: post
title: "The Story so Far: Part III, Assemblies and Alignments [WIP]"
---

<p class="message"><b>This post also isn't finished yet. Consider coming back later?</b></p>

At the end of my [previous post]({{ page.previous.url }}), I was left in a crumpled heap after learning the hard way that
not only is the `FASTQ` format specification somewhat non-existent but that I am also poor at reading reports I myself
generated; *wasting*[^1] time attempting to solve a problem that simply did not exist.

So, with our two large FASTQ files containing ~195 million pairs of reads sequenced from the contents of ~60 extracted
limpet gastroinestinal tracts, now *stringently* confirmed to have the same number of sequences. I could finally begin to
look at the next station on my whistle-stop tour of the average metagenomics pipeline.

# Assembly
Consider a completed jigsaw, depicting whatever you like. I'm going to select Aberystwyth University's
Director of Postgraduate Study, but you can choose whatever.

![]({{ site.url }}/public/posts/so-far-p3/rrz_jigsaw_lg.jpg)

Simplistically; single **reads** can be thought of as puzzle pieces which share some overlap with other
locally-similar puzzle pieces. We use these similarities to group neighbouring pieces and **assemble**
them in to longer runs or clusters of pieces (called **contigs**) which in turn share some overlap with other
runs or clusters of pieces that we have constructed elsewhere on the table, which can also be assembled.
Hopefully, we have enough information (similarity and overlaps) as well as all of the pieces to complete
the jigsaw that is our genome of interest[^2].

If we had little (or no) information about similarity (*e.g.* no images on the pieces, *i.e.* genetically
dissimilar reads) or overlap (*e.g.* the pieces are merely tiles with no tabs or blanks[^3], *i.e.* reads
which share no overlap) you can imagine how things become somewhat more difficult. Especially so if pieces
are missing.

## velvet
### BIGASSEMBLY
### VBIGASSEMBLY

* * *

# tl;dr

[^1]: I use this word lightly as it was something to learn nonetheless.

[^2]: Though in reality we isolate and attempt to assemble all the 'edge' pieces first, before spending hours staring at a picture of what the finished puzzle looks like in the first place, scrabbling to organise pieces in some form of order and the whole genomes-as-a-jigsaw metaphor quickly unravels itself.

[^3]: [I had no idea either.](http://english.stackexchange.com/questions/47667/what-do-you-call-the-interconnecting-bits-of-a-puzzle-piece-in-english)
