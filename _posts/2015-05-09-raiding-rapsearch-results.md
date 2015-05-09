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
does the top-ranking hit represent a hydrolase? How is "top-ranking" defined? Let's read the manual[^3]:

<blockquote>[...] preference is given to the db quality first, than the bit score and finally the lenght of annotation, the one with the highest values is kept
<footer><a href="http://pythonhosted.org/mgkit/scripts/filter-gff.html#overlap-filtering">http://pythonhosted.org/mgkit/scripts/filter-gff.html#overlap-filtering</a></footer></blockquote>

Hmm, it sounds as though database quality is considered paramount by the default ranking device, ensuring
that SwissProt results take precedence over those from TrEMBL. However what if the SwissProt hit actually has
a lower bitscore? To check, I'll modify [the source](https://bitbucket.org/setsuna80/mgkit/src/f3b5aa7e65d1cc8870743a9c7492ccb2528b8417/mgkit/filter/gff.py?#cl-45)[^4] and construct a trivial example below[^5].

```python
# lambda a1, a2: min(a1, a2, key=lambda el: (el.dbq, el.bitscore, len(el)))
discard = lambda a1, a2: min(a1, a2, key=lambda el: (el[0], el[1], el[2]))
crap_sp_hit = [10, 39, 100]
good_tr_hit = [8, 60, 100]
discard(crap_sp_hit, good_tr_hit)
> [8, 60, 100]
```

Welp, the "better" TrEMBL originating annotation is discarded. It seems the default selection function
ensures database quality above all else. Ideally we'd like a metric that gives weight to sequences held
in SwissProt (to reflect their curation accuracy) but not so much that they are always chosen over better
hits from a "weaker" database. `filter-gff` does accept an optional `--choose-func` parameter which I will
now be investigating the behaviour of.

As an aside, it's interesting to note the very high number of hits for fungal-associated hydrolases,
especially after overlap filtering -- where for both SwissProt and TrEMBL, it has more
results than the bacteria and archaea. I wonder whether this is indicative of the host
contamination known to be in the data set, as I guess host associated sequences are more
likely to have hits in the fungal database as they are both at least eukaryotic.

Before moving on, I'll re-run the overlap filtering step with a less naive filtering function and report back.
I'm curious to see what other **overlap loss** is being incurred, are the overlapping hits similar in function
and taxonomy, or widely different?

To really consider whether or not the hits are representative of a hydrolase, we need to calculate
how much the hit covers the whole database target sequence. It's all well and good if we have a high
bitscoring hit to a hydrolase, but if it only covers a fraction of the whole sequence, that doesn't
necessarily bode well for a "real" hydrolase being on the contig. Unfortunately, the `m8`/`blast6`
output formats (such as that of `rapsearch`) do not give the length of the target sequences, electing
to only give the length of the hit region and its identity.

So the next step will be to index the `FASTA` files used to build the hydrolase databases, then for each hit
found by `rapsearch`: query the `FASTA` indices for the length of the target hydrolase and work out the proportion
of the target covered by the hit. I can then add these values to the `GFF` files (containing the hits)
and refer to them in my own hit discarding function when re-calling `filter-gff`. Easy.

I'm somewhat confused though, the proportion of a database sequence covered by a hit seems like a common
question to determine how "good" a hit really is? I spoke to Tom and he normally defines a good hit by its bitscore
and by the number of bases on the hit that actually matched the target sequence exactly
(`hit length * hit identity`), he tells me this is pretty common too.

I guess we want to ensure we're discovering "real" hydrolases by picking out annotations that cover
as much of a known hydrolase gene as possible. But it still seems strange to me that this sort of
metric isn't more commonly used, are we doing something special?

* * * 

# tl;dr
* We have data. Already it raises more questions than answers.
* MGKit's `filter-gff` annotation overlap discarder liberally discards annotations from a "weaker" database by default.
* You never quite get what you need from a program's output.

[^1]: This seems to be a fairly commonly used "cut-off" threshold when looking at hit quality.

[^2]: Using [filter-gff](http://pythonhosted.org/mgkit/scripts/filter-gff.html#overlap-filtering).

[^3]: This is a good reflex reaction.

[^4]: I did panic briefly when I saw the function minimises rather than maximises quality; only to read the documentation further and discover the function must return the annotation to be discarded, rather than selected. Phew.

[^5]: Our team uses a database quality score of `8` for **TrEMBL** and `10` for **Swissprot** annotations. Though the numbers don't particularly matter, just so long as `SP > TR`.
