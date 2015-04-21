---
layout: post
title: "The Story so Far: Part I, A Toy Dataset"
pubdir: "so-far-p1"
---

In this somewhat **long** and long overdue post; I'll attempt to explain the work done so far and an overview
of the many issues encountered along the way and an insight in to why doing science is much harder than
it ought to be.

<p class="message">This post got a little longer than anticipated, so I've sharded it like everything else
around here.</p>

* * *

## In the beginning...
To address my lack of experience in handling metagenomic data, I was given a small[^1] dataset to play with.
Immediately I had to adjust my interpretation of what constitutes a "small" file. Previously the largest
single input I've had to work with would probably have been the human reference genome ([GRCh37](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/)) which as a FASTA[^2] file
clocks in at around a little over 3GB[^3].

Thus imagine my dismay when I am directed to the directory of my input data and find 2x**42GB** file.  
Together, the files are 28x the size of the largest file I've ever worked with...

### So, what is it?
Size aside, **what** are we even looking at and how is there so much of it?

The files represent approximately 195 million read pairs from a nextgen[^4] sequencing run, with each file holding
each half of the pair in the FASTQ format. The dataset is from a previous IBERS PhD student and
was introduced in a 2014 paper titled [Metaphylogenomic and potential functionality of the limpet Patella pellucida's gastrointestinal tract microbiome \[Pubmed\]](http://www.ncbi.nlm.nih.gov/pubmed/25334059). According
to the paper over 100 Blue-rayed Limpets (*Patella pellucida*) were collected from the shore of Aberystwyth, placed
in to tanks to graze on Oarweed (*Laminaria digitata*) for one month. 60 were plated, anesthetized
and aseptically dissected; vortexing and homogenizing the extracted digestion tracts before repeated
filtering and final centrifugation to concentrate cells as a pellet. The pellets were then resuspended and
DNA was extracted with a soil kit to create an Illumina paired-end library.

The paper describes the post-sequencing data handling briefly: the net result of 398 million reads which
were quality processed using `fastq-mcf`; to remove adaptor sequences, reads with quality lower than 20 and reads shorter than 31bp. The first 15bp of each read were also truncated[^5]. It was noted the remaining 391 million reads were heavily contaminated with host-derived sequences and thus insufficient for functional analysis.

My job was to investigate to what extent the contamination had occurred and to investigate whether
any non-limpet reads could be salvaged for functional analysis.

Let's take a closer look at the format to see what we're dealing with.

#### FASTQ Format
FASTQ is another text based file format, similar to FASTA but also stores quality scores for
each nucleotide in a sequence[^6]. Headers are demarcated by the `@` character instead of `>`
and although not required tend to be formatted strings containing information pertaining to the
sequencing device that produced the read.
Sequence data is followed by a single `+` on a new line, before a string of quality scores (encoded as
ASCII characters within a certain range, depending on the quality schema used) follows on
another new line:

```
@BEEPBOOP-SEQUENCER:RUN:CELL:LANE:TILE:X:Y 1:FILTERED:CONTROL:INDEX
HELLOIAMASEQUENCEMAKEMEINTOAPROTEINPLEASES
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
```

This example uses Illumina 1.8 sequence identifiers and quality scores, the same as those found
in the dataset. The quality string represents increasing quality from 0 (worst) to 41 (best) left to right.
Taking the first read from the first file as an actual example, we get a better look at a "real-world" sequence header:

```
@HWI-D00173:21:D22FBACXX:5:1101:1806:1986 1:N:0:CTCTCT
TTGTGTCAAAACCGAACAACATGACAATCTTACTTGCCTGGCCCTCCGTCCTGCACTTCTGGCATGGGGAAACCACACTGGGGGC
+
IIIAEGIIIFIIIEGIFFIIIFIFIIEFIIIIEFIIEFGCDEFFFFABDDCCCCCBBBBBBBBBBBBBB?BBBB@B?BBBBBBB5
```

So how are these files so large[^7]? Given each read record takes four lines 
(assuming reads short enough to not be spead over multiple lines -- which they are not
in our case) and each file contains around 195 million reads, we're looking at 780 million lines. Per file.

The maximum sequence size was 86bp and each base takes one byte to store, as well as corresponding
per-base quality information:

```
86 * 2 * 195000000
> 33540000000B == 33.54GB
```

Allowing some arbitrary number of bytes for headers and the `+` characters:

```
((86 * 2) + 50 + 1) * 195000000
> 43485000000B == 43.49GB
```

It just adds up. To be exact, both input files span 781,860,356 lines each -- meaning around
781MB of storage is used merely for newlines alone! It takes `wc` 90 seconds to count lines,
these files aren't small at all!


### Quality Control and Trimming

Although already done (as described by the paper), it's good to get an idea of how to run and interpret
basic quality checks on the input data. I used [`FASTQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
which outputs some nice HTML reports that can be archived somewhere if you are nice and organised.

<p class="message"><b>Top Tip</b><br />
It's good to be nice and organised because in writing this blog post I've been able
to quickly retreive the FASTQC reports from October and realised I missed a glaring problem
as well as a metric that could have saved me from wasting time.</p>

For an input FASTQ file, FASTQC generate a summary metrics table.
I've joined the two tables generated for my datasets below.

| Measure            | Value (\_1)                         | Value (\_2)                         |
|--------------------|-------------------------------------|-------------------------------------|
| Filename           | A3limpetMetaCleaned\_1.fastq.trim   | A3limpetMetaCleaned\_2.fastq.trim   |
| File type          | Conventional base calls             | Conventional base calls             |
| Encoding           | Sanger / Illumina 1.9               | Sanger / Illumina 1.9               |
| Total Sequences    | 195465089                           | 195465089                           |
| Filtered Sequences | 0                                   | 0                                   |
| Sequence length    | 4-86                                | 16-86                               |
| %GC                | 37                                  | 37                                  |

Here I managed to miss two things:

* Both files store the same number of sequences (which is expected as the sequences are paired), something that I apparently forget about shortly...
* Both files do not contain sequences of uniform length, neither do the non-uniform lengths have the same range, meaning that some pairs will not align to correctly as they cannot overlap fully...

FASTQC also generates some nice graphs, of primary interest, per-base sequence quality over the length of a read.

\_1            |  \_2
:-------------------------:|:-------------------------:
![]({{ site.url }}/public/posts/{{ post.pubdir }}/pbq1.png)  |  ![]({{ site.url }}/public/posts/{{ post.pubdir }}/pbq2.png)

```bash
# Command: LC_ALL=C grep -c '^@' $FILE
196517718   /ibers/ernie/repository/genomics/mts11/READS/Limpet-Magda/A3limpetMetaCleaned_1.fastq.trim
196795722   /ibers/ernie/repository/genomics/mts11/READS/Limpet-Magda/A3limpetMetaCleaned_2.fastq.trim
```

Let's take another look at the valid range of characters for the Illumina 1.8+ quality scores:
```
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
```
**`@`**

...for some reason someone thought it would be a great idea to allow the `@` character to
be available for use in quality strings...
...I'd written a lazy parser that just looked for these instead of checking for the correct
formatting of sequencing `+` quality, following by a new sequence starting with the `@` symbol...

Unfortunately this led me to believe the files had been incorrectly trimmed (with some reads in Set A
not appearing in Set B) launching me head first in to my first large scale
computing problem. Given a set of 195,000,000 reads, and another set of mostly similar 195,000,000 reads
 Can we distinguish what reads appear in one set and not the other, given that, can we then create an output file
 that only contains the intersection of the two sets?
 
 At first I tried to use Python dictionaries to store all the read data seen in both of the input files. It look a while to load. A long while... Turns out that Python requires almost quarter of an hour to merely read through the files and count the number of lines.  I was thinking something was a little off about this, then considered that there are XX billion lines in the file.... Even grep takes 2 minutes.
 
 
```bash
# Command: LC_ALL=C grep -c '^@HWI-D' $FILE
195465089   /ibers/ernie/repository/genomics/mts11/READS/Limpet-Magda/A3limpetMetaCleaned_1.fastq.trim
195465089   /ibers/ernie/repository/genomics/mts11/READS/Limpet-Magda/A3limpetMetaCleaned_2.fastq.trim
```

Better! The number of sequence headers match and I'm happy to assume they are all correctly paired
as that would require something pretty screwy to be going on\*\*\*\*\*\*.
Although even this check is only *just about* robust, assuming that all reads are from the same suite
of machines addressed as `HWI-DXXXXX` (which I do know for sure) 
and secondly that `H`, `I` and `-` are all valid quality score
characters and so `W` is the only character preventing further accidental matches to quality strings.



* * *

There's two main issues of size here:
* Jobs that are large because the inputs are large,
* Jobs that are large because the inputs are tiny but there are thousands of them...

* * *
# tl;dr
* Don't try and count the number of sequences in a FASTQ file by counting `@` characters.
* Prepend `LC_ALL=C` to commands like `grep` and `awk` if you don't need to support non-ASCII character spaces.
* Processing **massive** files takes time (more than a minute) and there's nothing wrong with that.
* Practice reading
* Try to actually read quality reports, then read them again. Then grab a coffee and read them for a third time before you do anything.

* * *

[^1]: Now realised to be a complete misnomer, both in terms of size and effort.

[^2]: A text based file format where sequences are delimited by `>` and a sequence name [and|or] description,
    followed by any number of lines containing nucleotides or amino acids (or in reality, whatever you fancy):

    ```bash
    >Example Sequence Hoot Factor 9 | 00000001
    HELLOIAMASEQUENCE
    BEEPBOOPTRANSLATE
    MEINTOPROTEINS
    >Example Sequence Hoot Factor 9 | 00000002
    NNNNNNNNNNNNNNNNN
    ```  
    
    Typically sequence lines are of uniform length (under 80bp), though this is not a requirement of the format.
    The [NCBI](http://www.ncbi.nlm.nih.gov/) suggest formats for the header (single line descriptor,
    following the '>' character) though these are also not required to be syntactically valid.

[^3]: Stored as text we take a byte for each of the 3 billion nucleotides as well as each newline
    delimiter and an arbitrary number of bytes for each chromosome's single line header.

[^4]: Seriously, can we stop calling it nextgen yet?

[^5]: I'm unsure why, from a recent internal talk I was under the impression we'd normally trim the first "few" bases (3-5bp, maybe 8bp if there's a lot of poor quality nucleotides) to try and improve downstream analysis such as alignments (given the start and end of reads can often be quite poor and not align as well as they should) but 15bp seems excessive. It also appears the ends of the reads were not truncated.

[^6]: Which actually wouldn't be that much of a surprise.

[^7]: Or small, depending on whether you've adjusted your world view yet.
