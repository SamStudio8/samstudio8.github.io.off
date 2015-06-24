---
layout: post
title: "Sanger Sequel"
---

In a change to scheduled programming, days after touching down from my holiday (which needs a post of its own)
I moved[^4] to spend the next few weeks back at the
[Wellcome Trust Sanger Institute](https://www.sanger.ac.uk/) in Cambridgeshire. I interned here
previously in 2012 and it's still like working at a science-orientated Google thanks to the
overwhelming amount of work being done and the crippling inferiority complex that comes
from being surrounded by internationally renowned scientists. Though at least I'm not acquiring
significant mass from free food.

My aim here is two fold and outlined below. Though of course, I'm still on the books back at Aberystwyth and it
would be both naughty and cruel of me to leave my [newly acquired data]({% post_url 2015-06-01-rapsearch-returns %})
cold, alone and untouched until I get back.

## Aims
### (1) Produce a Sequel
My [undergraduate dissertation](https://github.com/SamStudio8/frontier-dissertation) was titled:
*Application of Machine Learning Techniques to Next Generation Sequencing Quality Control* and
worked in collaboration with some colleagues from my previous placement at the Sanger Institute[^1].
The project was to build a machine learning framework capable of improving detections of "bad" samples
by first characterising what it meant to be a bad sample.

In short, the idea was to repeatedly push a large number of samples (each known to have individually
passed or failed some internal quality control mechanism) through some analysis pipeline, holding out
a single sample out from the analysis in turn. The difference to a known result would then be
calculated and samples would be re-classified as good or bad based on whether the accuracy of a
particular run was increase or decreased in their absence.

Ultimately the scope was too large and the tools too fragile to complete the end-goal in the time
that I had (though it still achieved 90% and won an award, so one can't complain too much) but
we still have the data and while I am here it would be interesting to try and pick up where we left off.
I expect to do battle with the following tasks over the next few days:

**Recall in detail what we were doing and figure out how far we got**  
*i.e.* Dig out the thesis, draw some diagrams and run `ls` everywhere.

**Confirm the `Goldilocks` region**  
Due in part to the short time that I had to complete this project the first time around -- a constraint
I still have -- I authored a tool named [`Goldilocks`](https://goldilocks.readthedocs.org) to "narrow down"
my analysis from a whole genome to just a 1Mbp window. It would be worth ensuring the latest version of
`Goldilocks` (which has long fixed [some bugs I would really like to forget](https://github.com/SamStudio8/goldilocks/commit/b6b1f6f560202d6e33df3bfcec1d48a35fe8c6c0)) still returns the same results as it did when I
was doing my thesis.

**Confirm data integrity**  
The data has been sat around in the way of actual in-progress science for the best part of a year
and has possibly been moved or "tidied away". It would be worth ensuring all the files are actually
intact and for the sake of completeness revisit how those files came to be and regenerate them. This will
encompass ensuring the `Goldilocks` region for each sample was correctly extracted. I recall the samples
were made up of two studies and we may have decided not to pursue one of them due to differences in sequencing[^2].
I also recall having some major trouble with needing to re-align the failed samples to a different reference: these
samples having failed, were not subjected to all the processing of their QC approved counterparts, which we'll need
to apply ourselves manually, presumably painstakingly.

**Prepare data for the pipeline**  
The nail in the coffin for the first stab at this project was the data preparation: `samtools merge` was just
woefully slow in handling the scale of data that I had, in particular struggling to merge many thousands of files
at once. A significant amount of project time was spent tracking and patching memory leaks and contributing other
functionality (more on this in a moment) that left me with little time at the end to actually push the data
through the pipeline and get results. `samtools` has undergone some rapid improvements since and I suspect
this step will no longer pose such a hurdle.


### (2) Contribute to `samtools`
As I briefly alluded, during the course of my undergraduate dissertation
I authored several pull requests to a popular open-source bioinformatics toolkit known as `samtools`,
which was initially created and continues to be maintained right here at the Sanger Institute.
In particular, these pull requests improved documentation and patched some memory leaks for `samtools merge`
and also added naive header parsing for input file metadata to be organised into basic structures for
much more efficient iterative access later; significantly improving the time performance of `samtools merge`.

Header parsing has been a long sought after feature for `samtools` but none of the core
maintainers had the time to put aside to take a good look at the RFC I had submitted. Now I'm in-house
and I put a face to a username, catching the most recent `samtools` steering meeting off-guard, I've been
tasked to try and get this done before I leave at the end of July.

No pressure.

* * *

#tl;dr
* I live in Cambridge until the end of July, please don't try and find me in my office.[^3]


[^1]: As soon as I had my foot in the door I refused to take it away.

[^2]: The sequencing was conducted at different depths between the two standalone studies and it was suspected this may introduce some bias I didn't want to deal with.

[^3]: Not that I'm ever in there, anyway.

[^4]: For what is *at least* the 10th major house move I've made since heading out to university.
