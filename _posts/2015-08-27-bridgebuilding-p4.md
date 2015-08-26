---
layout: post
title: "The Tolls of Bridge Building: Part IV, Mysterious Malformations [WIP]"
---

Following a short hiatus on the sample *un*improvement job which may or may not
have been halted by `vr-pipe` inadvertently [knocking over a storage node]({{ page.previous.url }})
at the Sanger Institute, our jobs burst back in to life
only to fall at the final hurdle of the first pipeline of the `vr-pipe`
workflow. Despite my lack of tweed hat and pipe, it was time to play
bioinformatics Sherlock Holmes.

Without `mercury` access to review reams of logs myself, the `vr-pipe` web interface
provides a single sample of output from a random job in the state of interest.
Taking its offering, it seemed I was about to get extremely familiar with lanelet `7293_3#8`.
The error was at least straightforward...

### Accounting Irregularities

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

Though, my first instinct was to check whether the numbers stated in the message were
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

Every single read in the "missing" `285,720` can be accounted for as filtered by this `MalformedReadFilter`.
This definitely isn't expected behaviour, as I said already, we have ways of marking reads as QC failed, supplementary or unmapped without just discarding them wholesale. Our pipeline is not supposed to drop data.

* * *

#tl;dr
* `GATK` has a `MalformedReadFilter`
