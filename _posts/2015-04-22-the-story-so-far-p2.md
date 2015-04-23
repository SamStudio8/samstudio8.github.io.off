---
layout: post
title: "The Story so Far: Part II, An Unnecessary Tangent [WIP]"
---

<p class="message"><b>This post isn't finished yet. Consider coming back later?</b></p>

## An Unnecessary Tangent
<p class="message">Whilst a learning exercise, I feel it necessary to re-iterate that this tangent could easily have been avoided had I only paid a little more attention to the `FASTQC` reports...<br /><br />
If you are after a project overview, you can skip this interluding derailment and head to the next part.</p>

### Hello Python
I jumped to write a Python script that reads in each file by turn, adds each sequence header to a dictionary and stores
the byte-location of that read in both files (from `file_handle.tell()`) in a two element list. Many hours after
executing this script I was confused that I was still waiting, progress was happening but becoming increasingly
sluggish. I'd saturated the dictionary. ~196 million keys is a lot.

I started to wonder about things I'd not noticed or taken in to account before -- how long does it even take to just
handle the file? Iterating over one of the files in Python with `readline()` took just shy of quarter of an hour[^1].
Wat?

```python
while line:
    if line[0] == '@':
        count += 1
    line = fastq_1_fh.readline()
print count 
```

This was unsettling. This is the sort of thing I've not needed to worry about before. Files just get handled and
then things happen? I've never experienced such a **wait**! I scrabbled for a calculator, to take 15 minutes to
handle ~781 million lines, we must be doing around 867,777 lines every second, that's not bad.
Suddenly the scale of what I was trying to do set in: it's not that Python is *slow* (although neither is it the
fastest method), there really is just a **significant** amount of data here. Over three quarters of a billion
lines; that's twelve lines (or rather, three whole FASTQ records) for every person in the UK[^2].

This set the baseline for improvement, we couldn't possibly do what I wanted to do in less than 30 minutes as file
handling alone took this much time.

### Oddities with `grep` and `awk`

### MySQL

### Mootness
...
```bash
LC_ALL=C grep -c '^@HWI-D' $FILE
```
```
195465089   A3limpetMetaCleaned_1.fastq.trim
195465089   A3limpetMetaCleaned_2.fastq.trim
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
* Prepend `LC_ALL=C` to commands like `grep` and `awk` if you don't need to support non-ASCII character spaces.
* Python loads file in blocks to speed up file handling so unless you control the iterator yourself (manually calling `readline()`, `tell()` will respond with the location of the end of the current block, which is a bit useless if you are trying to build an index of where things are...
* Processing **massive** files takes time (more than a minute) and there's nothing wrong with that.

* * *

[^1]: In reality I wouldn't go so far as describing this as a "benchmark", my methodology is a little flawed. I was using the login/head node (which experiences varying levels of load) and runs Python 2.6 by default. I also made no real attempt to flush cache lines. Despite this, I demonstrated reproducible behaviour and thought it worth of mention.

[^2]: Or 65 records for each person here in Wales.
