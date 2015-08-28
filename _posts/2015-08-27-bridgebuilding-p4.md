---
layout: post
title: "The Tolls of Bridge Building: Part IV, Mysterious Malformations"
---

Following a short hiatus on the sample *un*-improvement job which may or may not
have been halted by `vr-pipe` inadvertently [knocking over a storage node]({% post_url 2015-07-31-bridgebuilding-p3 %})
at the Sanger Institute, our **837** *non-33* jobs burst back in to life
only to fall at the final hurdle of the first pipeline of the `vr-pipe`
workflow. Despite my lack of deerstalker and pipe, it was time to play
bioinformatics Sherlock Holmes.

Without `mercury` access to review reams of logs myself, the `vr-pipe` web interface
provides a single sample of output from a random job in the state of interest.
Taking its offering, it seemed I was about to get extremely familiar with lanelet `7293_3#8`.
The error was at least straightforward...

### Accounting Irregularities[^0]

```
... [step] failed because 48619252 reads were generated in the output bam
file, yet there were 48904972 reads in the original bam file at [...]
/VRPipe/Steps/bam_realignment_around_known_indels.pm line 167
```

The first thing to note is that the error is raised in `vr-pipe` itself,
not in the software used to perform the step in question -- which happened to be `GATK`,
for those interested. `vr-pipe` is open source software,
[hosted on Github](https://github.com/VertebrateResequencing/vr-pipe). The deployed release of
`vr-pipe` used by the team [is a fork](https://github.com/wtsi-hgi/vr-pipe) and so the
[source raising the error](https://github.com/wtsi-hgi/vr-pipe/blob/hgi-release/modules/VRPipe/Steps/bam_realignment_around_known_indels.pm#L166)
is available to anyone with the stomach to read Perl:

```perl
my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
my $actual_reads = $out_file->num_records;
        
if ($actual_reads == $expected_reads) {
    return 1;
}
else {
    $out_file->unlink;
    $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output bam file, yet there were $expected_reads reads in the original bam file");
}
```

The code is simple, `vr-pipe` just enforces a check that all the reads from the
input file make it through to the output file. Seems sensible. This is very much
desired behaviour as we never drop reads from BAM files, there's enough flags and
scores [in the SAM spec](https://samtools.github.io/hts-specs/SAMv1.pdf) to put reads
aside for one reason or another without actually removing them.
So we know that the job proper completed successfully (there was an output file sans
errors from `GATK` itself). So the question now is, *where did those reads go?*

Though before launching a full scale investiation, my first instinct was to question the error
itself and check the numbers stated in the message were
even correct. It was easy to confirm the original number of reads: `48,904,972` by checking
both the `vr-pipe` metadata I generated to submit the bridged BAMs to `vr-pipe`. I ran
`samtools view -c` on the bridged BAM again to be sure.

But `vr-pipe`'s standard over-reaction to job failure is nuking all the output files, so I'd
need to run the step myself manually on lanelet `7293_3#8`. A few hours later, `samtools view -c`
and `samtools stats` confirmed the output file really did contain `48,619,252` reads, a shortfall of `285,720`.

I asked Irina, one of our `vr-pipe` sorceresses with `mercury` access to dig out the
whole log for our spotlight lanelet. Two stub `INFO` lines sat right at the tail of the `GATK`
output shed immediate light on the situation...

### Expunged Entries

> INFO 18:31:28,970 MicroScheduler - 285720 reads were filtered out during the traversal out of approximately 48904972 total reads (0.58%)
> 
> INFO 18:31:28,971 MicroScheduler - -> 285720 reads (0.58% of total) failing MalformedReadFilter

Every single of the `285,720` "missing" reads can be accounted for by this `MalformedReadFilter`.
This definitely isn't expected behaviour, as I said already, we have ways of marking reads as QC failed, supplementary or unmapped without just discarding them wholesale. Our pipeline is not supposed to drop data.
The [`MalformedReadFilter` documentation](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_engine_filters_MalformedReadFilter.php) on the `GATK` website states:

> Note that the MalformedRead filter itself does not need to be specified in the command line because it is set automatically.

This at least explains the "unexpected" nature of the filter, but surely `vr-pipe` encounters files with bad
reads that need filtering all the time? My project can't be the first...
I figured I must have missed a pre-processing step, I asked around: "Was I supposed to do my own filtering?". I compared
the `@PG` lines documenting programs applied to *the 33* to see whether they had undergone different treatment,
but I couldn't see anything related to filtering, quality or otherwise.

I escalated the problem to Martin who replied with a spy codephrase:

> the MalformedReadFilter is a red herring

Both Martin and Irina had seen *similar* errors that were usually indicative of the wrong version of
`vr-pipe` [and|or] `samtools` -- causing a subset of flagged reads to be counted rather than all.
But I explained the versions were correct and I'd manually confirmed the reads as actually missing
by running the command myself, we were all a little stumped.

I read the manual for the filter once more and realised we'd ignored the gravity of the word "malformed":

> This filter is applied automatically by all GATK tools in order to protect them from crashing on reads that are grossly malformed.

We're not taking about bad quality reads, we're talking about reads that are incorrect in such a way
that it may cause an analysis to terminate early. My heart sank, I had a hunch. I ran the output from
lanelet `7293_3#8`'s `bridgebuilder` adventure through [`GATK PrintReads`](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php), an arbitrary tool
that I knew also applied the `MalformedReadFilter`. Those very same `INFO` lines were printed.

My hunch was right, the reads had been malformed all along.

### Botched Bugfix
I had a pretty good idea as to what had happened but I ran the inputs to `brunel`
(the final step of `bridgebuilder`) through `PrintReads` as a sanity check. This proved
a little more difficult to orchestrate than one might have expected, I had to falsify headers
and add those pesky missing `RG` tags that plagued us before.

The inputs were fine, as suspected, `brunel` was the malformer. My hunch? My [quick hack](https://github.com/SamStudio8/bridgebuilder/commit/e5a83b46d64f888c6dc7780779a2786d7849328e)
to [fix my last `brunel` problem]({% post_url 2015-07-31-bridgebuilding-p3 %}#translations) had come back
to bite me in the ass and caused an even more subtle `brunel` problem.

Indeed, despite stringently checking the bridged BAMs with
[five different tools]({% post_url 2015-07-23-bridgebuilding-p2 %}#finalchecks), successfully generating
an index and even processing the file with `Picard` to mark duplicate reads and `GATK` to detect and re-align
around indels, these malformed reads *still* flew under the radar -- only to be caught by a few
lines of Perl that a minor patch to `vr-pipe` happened to put in the way.

Recall that my `brunel` fix initialises the translation array with `-1`:

> **Initialise `trans[i] = -1`**  
> The only quick-fix grade solution that works, causes any read on a TID that has no translation to be regarded as "unmapped". Its TID will be set to "*" and the read is placed at the end of the result file. The output file is however, valid and indexable.

This avoided an awful bug where `brunel` would assign reads to chromosome 1 if their actual
chromosome did not have a translation listed in the user-input translation file.
In practice, the fix worked. Reads appeared "unmapped", their TID was an asterisk and they were listed at
the end of the file. The output was viewable, indexable and usable with `Picard` and `GATK`,
but *technically* **not** valid afterall.

To explain why, let's feed everybody's favourite bridged BAM `7293_3#8` to [`Picard ValidateSamFile`](http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile), a
handy tool that does what it says on the tin. `ValidateSamFile` linearly passes over each record of a
BAM and prints any associated validation warnings or errors to a log file[^1]. As the tool moved over the target,
an occassional warning that could safely be ignored was written to my terminal[^2]. The file seemed valid and
as the progress bar indicated we were nearly out of chromosomes, I braced.

As `ValidateSamFile` attempted to validate the unmapped reads, a firehose of errors (I was expecting
a lawn sprinker) spewed up the terminal and despite my best efforts I couldn't `Ctrl-C` to point away
from face.

### False Friends
There was a clear pattern to the log file. Each read that I had translated to unmapped triggered
four errors. Taking just one of thousands of such quads:

> ERROR: [...] Mate Alignment start should be 0 because reference name = *.  
> ERROR: [...] Mapped mate should have mate reference name  
> ERROR: [...] Mapped read should have valid reference name  
> ERROR: [...] Alignment start should be 0 because reference name = *.

Now, a bitfield of flags on each read is just one of the methods employed by the [BAM format](http://samtools.github.io/hts-specs/SAMv1.pdf) to perform validation
and filtering operations quickly and efficiently. One of these bits: `0x4`, indicates a read is
unmapped. Herein lies the problem, although I had instructed `brunel` to translate the TID (which
refers to the i-th `@SQ` line in the header) to `-1` (*i.e.* no sequence), I did not set the unmapped flag.
This is invalid (`Mapped read should have valid reference name`), as the read will appear as aligned,
but to an unknown reference sequence.

Sigh. A quick hack combined with my ignorance of the underlying BAM format specification was at fault[^4].

The fix would be frustratingly easy, a one-liner in `brunel` to raise the `UNMAPPED` flag (and to be
good, a one-liner to unset the `PROPER_PAIR` flag[^3]) for the appropriate reads. Of course, expecting
an easy fix jinxed the situation and the `MalformedReadFilter` cull persisted, despite my new
semaphore knowledge.

For each offending read, `ValidateSamFile` produced a triplet of errors:

> ERROR: [...] Mate Alignment start should be 0 because reference name = *.  
> ERROR: [...] MAPQ should be 0 for unmapped read.  
> ERROR: [...] Alignment start should be 0 because reference name = *.

The errors at least seem to indicate that I'd set the `UNMAPPED` flag correctly.
Confusingly, the [format spec](http://samtools.github.io/hts-specs/SAMv1.pdf) has the following point (parentheses and emphasis mine):

> Bit 0x4 (**unmapped**) is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, and bits 0x2 (**properly aligned**), 0x100 (secondary), and 0x800 (supplemental).  

It would seem that canonically, the alignment start (`POS`) and mapping quality (`MAPQ`) are untrustworthy
on reads where the unmapped flag is set. Yet this triplet appears exactly the same number of times
as the number of reads hard filtered by the `MalformedReadFilter`.

I don't understand how these reads could be regarded as "grossly malformed" if even the format specification
acknowledges the possibility of these fields containing misleading information. Invalid yes, but grossly malformed? No.
I just have to assume there's a reason `GATK` is being especially anal about such reads, perhaps developers simply (and rather fairly)
don't want to deal with the scenario of not knowing what to do with reads where the values of half the fields can't be trusted[^6].
I updated `brunel` to set the alignment positions for the read and its mate alignment to position 0 in retaliation.

I'd reduced the `ValidateSamFile` erorr quads to triplets and now, the triplets to duos:

> ERROR: [...] Mate Alignment start should be 0 because reference name = *.  
> ERROR: [...] Alignment start should be 0 because reference name = *.

The alignment positions are non-zero? But I just said I'd fixed that? What gives?

I used `samtools view` and grabbed the tail of the BAM, those positions were indeed non-zero.
I giggled at the off-by-one error, immediately knowing what I'd done wrong.

The BAM spec describes the alignment position field as:

> **POS**: 1-based leftmost mapping POSition of the first matching base.

But under the hood, the `pos` field is
[defined in `htslib`](https://github.com/samtools/htslib/blob/f79045536b9bb89666bedf7ee665503a0d9bad2a/htslib/sam.h#L138)
as a 0-based co-ordinate, because that's how computers work. The 0-based indices are then converted by just adding 1
whenever necessary. Thus to obtain a value of 0, I'd need to set the positions to -1. This off-by-one mismatch
is a constant source of embarrassing mistakes in bioinformatics, where typically 1-based genomic indices must be reconciled with 0-indexed data structures[^5].

Updating `brunel` once more, I ran the orchestrating `Makefile` to generate what I really hope to be the final
set of bridged BAMs, ran them through Martin's `addreplacerg` subcommand to fill in those missing `RG` tags
and then each went up against each of the now-six final check tools (yes, they validated with `ValidateSamFile`).
I checked index generation, re-ran a handful of sanity checks (mainly ensuring we hadn't lost any reads), re-generated
the manifest file before finally asking Irina for what I really hope to be the last time, to reset my `vr-pipe` setup.

I should add, Martin suggested that I use `samtools fixmate`, a subcommand designed for this very purpose.
The problem is, for `fixmate` to know where the mates are, they must be adjacent to each-other; that is, sorted
by name and not by co-ordinate. It was thus cheaper both in time and computation to just re-run `brunel` than to
name-sort, `fixmate` and coord-sort and current bridged BAMs.

<a name="victory"></a>
### Victorious `vr-pipe` *Update: Hours later*
Refreshing the `vr-pipe` interface [intermittently](https://www.youtube.com/watch?v=8UhONY3-1os), I was waiting
to see a number of inputs larger than 33 make it to the end of the first pipeline of the workflow: represented
by a friendly looking green progress bar.

Although the number currently stands at 1, I remembered that all of *the 33* had the same
leading digit and that for each progress bar, `vr-pipe` will offer up metadata on one sample job. I hopefully
clicked the green progress bar and inspected the metadata, the input lanelet was not in the set of the 33
already re-mapped lanelets.

I'd done it. After almost a year, I've led my lanelet Lemmings to the end of the first level.

* * *

#tl;dr
* `GATK` has an anal, automatated and aggressive `MalformedReadFilter`
* `Picard ValidateSamFile` is a useful sanity check for SAM and BAM files
* Your cleverly simple bug fix has probably swapped a critical and obvious problem for a critically unnoticeable one
* Software is not **obliged** to adhere to any or all of a format specification
* Off by one errors continue to produce laughably simple mistakes
* I should probably learn the BAM format spec inside out
* Perl is still pretty grim

[^0]: Those reads were just resting in my account!

[^1]: Had I known of the tool sooner, I would have employed it as part of the extensive `bridgebuilder` quality control suite.

[^2]: Interestingly, despite [causing jobs to terminate]({% post_url 2015-07-31-bridgebuilding-p3 %}#vanish-rg), a read missing an `RG` tag is a warning, not an error.

[^3]: Although not strictly necessary as per the specification (see below), I figured it was cleaner.
    
   > Bit 0x4 (**unmapped**) is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, and bits 0x2 (**properly aligned**), 0x100 (secondary), and 0x800 (supplemental).  
   
   Parentheses and emphasis mine.

[^4]: Which only goes to further my *don't trust anyone, even yourself* mantra.

[^5]: When I originally authored [`Goldilocks`](https://goldilocks.readthedocs.org), I tried to be clever and make things easier for myself, electing to use a 1-based index strategy throughout. This was partially inspired by FORTRAN, which features 1-indexed arrays. In the end, the strategy caused more problems than it solved and I had to carefully had to tuck my tail between my legs and [return to a 0-based index](https://github.com/SamStudio8/goldilocks/commit/6915d8f0347d7e71f947d84f1f381816b0d8a1e7).

[^6]: Though, [looking at the documentation](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_engine_filters_MalformedReadFilter.php), I'm not even sure *what* aspect of the read is triggering the hard filter. The three additional command line options described don't seem to be related to any of the errors raised by `ValidateSamFile` and there is no explicit description of what is considered to be "malformed".
