---
layout: post
title: "Raiding `rapsearch` Results"
---

Finally. After all the trouble I've had [trying to scale `BLAST`]({% post_url 2015-04-27-what-am-i-doing %}),
running [out of disk space]({% post_url 2015-04-25-scratch %}),
[database accounting irregularities]({% post_url 2015-04-24-trembling %})
and [investigating `boost::archive::archive_exception`]({% post_url 2015-05-04-aligned-annihilation-sequel %}) when using `rapsearch`, **we have data**.

Thanks to the incredible speed of `rapsearch`, what I've been trying to accomplish over the past few months
with `BLAST` has been done in mere hours without the hassle of database or contig sharding. Quantifying the
accuracy of `rapsearch` is still something our team is working on, but Tom's initial results suggest comparable
performance to `BLAST`. For the time being at least, it means I can get things done.

As [previously described]({% post_url 2015-05-01-aligned-annihilation %}), I extracted bacterial, archaeal and
fungal associated hydrolases from both the SwissProt (manually curated)
and TrEMBL (automatically annotated) databases. The tables below summarise the number of "raw" hits from `rapsearch`,
the number of hits remaining after discarding hits with a bitscore of less than 40[^1], followed by the hits remaining
after selecting for the "best" hit for cases where hits overlap by 100bp[^2] or more.

## SwissProt
| Taxa         | Raw SP  | Bitscore Filter | Overlap Filter |
|--------------|---------|-----------------|----------------|
|Bacteria [2]  | 127,403 | 58,715          | 1,649          |
|Archaea [2157]| 17,141  | 2,326           | 341            |
|Fungi [4751]  | 78,083  | 34,180          | 3,222          |
|**Total**     | 222,627 | 95,221 (42.77%) | 5,212 (5.47%, *Raw*:2.34%)|

## TrEMBL
| Taxa         | Raw TR  | Bitscore Filter | Overlap Filter |
|--------------|---------|-----------------|----------------|
|Bacteria [2]  | 683,307 | 392,791         | 6,810          |
|Archaea [2157]| 79,738  | 34,950          | 1,486          |
|Fungi [4751]  | 345,160 | 190,936         | 7,379          |
|**Total**     | 1,108,205 | 618,677 (55.83%)| 15,675 (2.53%, Raw:1.41%)|


## Merged

| Taxa         | Raw All  | Bitscore Filter | Overlap Filter |
|--------------|----------|-----------------|----------------|
|**All**       | 1,330,832| 713,898 (53.64%)| 12,194 (1.71%, Raw: 0.92%)|

Initially, I had merged all the hits from both databases and all three taxa together to create
a super-hit list; yielding just over 12k reasonable quality (`>=BQ40`) hits to
play with after both filtering steps. However, I became concerned with what I coin
**overlap loss**: a significant number of hits were discarded in overlapping regions.
94.53% and 97.47% of our bitscore filtered hits were lost by overlap for
SwissProt and TrEMBL lines respectively!

I suspect due to the condensed nature of the assembly (around 16.7 billion base pairs from the
raw reads aligning to an assembly consisting of just 433 million base pairs, an average coverage of
~38x) there is likely to be a lot of overlap occurring in 100bp windows. The question is, how well
does the top-ranking hit represent a hydrolase? How is "top-ranking" defined?

<blockquote>[...] preference is given to the db quality first, than the bit score and finally the lenght of annotation, the one with the highest values is kept
<footer><a href="http://pythonhosted.org/mgkit/scripts/filter-gff.html#overlap-filtering">http://pythonhosted.org/mgkit/scripts/filter-gff.html#overlap-filtering</a></footer></blockquote>

(selected by database quality first [SwissProt takes precedence over TrEMBL)

With this in mind, rather than continuing
the pipeline with the 12k results, I'll stop here to further investigate what overlapping results
are being discarded first.

It's interesting to note the very high number of hits for fungal-associated hydrolases,
especially after overlap filtering -- where for both SwissProt and TrEMBL, it has more
results than the bacteria and archaea. I wonder whether this is indicative of the host
contamination known to be in the data set, as I guess host associated sequences are more
likely to have hits in the fungal database as they are both at least eukaryotic.

* * * 

# tl;dr
* We have data. Already it raises more questions than answers.

[^1]: This seems to be a fairly commonly used "cut-off" threshold when looking at hit quality.

[^2]: Using [filter-gff](http://pythonhosted.org/mgkit/scripts/filter-gff.html#overlap-filtering).
