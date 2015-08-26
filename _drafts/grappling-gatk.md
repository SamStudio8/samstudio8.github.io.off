## Claw
## Extensions
> Invalid command line: The GATK reads argument (-I, --input_file) supports only BAM/CRAM files with the .bam/.cram extension

## SequenceDictionary

...Claw (0)
## ReorderSam
## Permissions
> Unable to create a temporary BAM schedule file.  Please make sure Java can write to the default temp directory or use -Djava.io.tmpdir= to instruct it to use a different temp directory instead.

## Claw (TID)
## Parameters to GenomeLocParser are incorrect

## MAX_FILE_HANDLES_FOR_READ_ENDS_MAP
> Exception in thread "main" htsjdk.samtools.SAMException: [...].tmp not found  
> [...]  
> Caused by: java.io.FileNotFoundException: [...].tmp (Too many open files)

`MAX_FILE_HANDLES_FOR_READ_ENDS_MAP` must be <= `ulimit -n`
```
