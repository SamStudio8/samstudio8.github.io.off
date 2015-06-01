---
layout: post
title: "`rapsearch` Returns"
---

Following completion of my [most recent side-quest]({% post_url 2015-05-18-playing-phylogenetic-hide-seek %}) to find
a little more about who the protozoa actually are and where they live in the context of UniProt, I now had a starting
point to append to my archive of hydrolase records. I had already shown that around
[1,500 Ciliophora-associated hydrolases](http://www.uniprot.org/uniprot/?query=taxonomy%3A%22Ciliophora+%5B5878%5D%22+AND+ec%3A3.*&sort=score)
could be extracted from UniProt, but before continuing my hunt for relevant protozoans, I wanted to run a quick sanity check.

As of May 2015, querying for *all* records in UniProt (*i.e.* either the manually-curated **SwissProt** or automatically-annotated **TrEMBL**)
which are assigned an EC (Enzyme Classification) of `3.*` yields just over [1.4M results](http://www.uniprot.org/uniprot/?query=ec%3A3.*&sort=score).
Yet [in a recent introduction to my newly extracted databases]({% post_url 2015-05-01-aligned-annihilation %}),
my bacterial-associated hydrolases totalled (`reviewed + unreviewed`) 1.1M records -- so we've pulled out
almost 80% of all the hydrolases in UniProt already!

[I mentioned previously]({% post_url 2015-04-27-what-am-i-doing %}) that we had started using `rapsearch` as an
alternative to `BLAST` due to its execution speed. In fact `rapsearch` was capable of searching through all ~700K
limpet contig sequences (totalling ~433 megabases) against the largest of the hydrolase databases I had created
(~1.1M sequences, 408 megabases) without the need to shard either the contigs, or the database --
in a matter of **hours** as opposed to *weeks* (or as it has felt, forever).

Given the primary reason for taking particular taxa-associated records was to reduce cluster time, but
clearly we're able to adequately process the vast majority of records in a more than reasonable time.
Thus there doesn't appear to be a reason against just making a superdatabase of all the hydrolases in UniProt,
a "few" extra sequences seems somewhat moot in terms of computational complexity...

Thus I tabulate below the results of executing `rapsearch` over the limpet contigs for all hydrolases in both
SwissProt and TrEMBL:

| Database Source | #Records | #Nucleotides   |Execution Time (Max. GB RAM) | Raw Hits | Bitscore Filter | Overlap Filter |
|-----------------|----------|----------------|-----------------------------|----------|-----------------|----------------|
|SwissProt        |64,521    |26,516,938      |0:35:03 (95.37)              |604,867   |224,394 (37.10%) |13,706 (6.11%, *Raw*:2.27%) |
|TrEMBL           |1,335,692 |545,226,198     |3:39:30 (129.32)             |1,975,908 |979,220 (49.56%) |33,599 (3.43%, *Raw*:1.70%) |
|**Total**        |1,400,213 |571,743,136     |4:14:33 (224.69)             |2,580,775 |1,203,614 (46.64%)|35,756 (2.98%, *Raw*:1.39%)|

