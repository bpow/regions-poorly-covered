# regions-poorly-covered

This is a wrapper around the GenomeAnalysisTK (GATK) DepthOfCoverage (DoC) walker to answer the question:

> Which regions from a large-scale genetic sequencing experiment are _consistently_ covered poorly?

Or, phrased more specifically:

> Given a set of ".bam" files and set of intervals of interest, which regions
> have "callable coverage" (filtering for mapping and base quality values) of
> less than 20x in at least 90% of the samples?

## Installation/Requirements

1. python version 2.7 (this would be trivial to re-write for an earlier python version-- would
   just need to replace argparse with optparse; or for a python3...)
2. A java version supported by GATK
2. `GenomeAnalysisTK.jar`: This is not redistributable, so you would have to register and accept
   the license agreement at [the GATK webpage](https://www.broadinstitute.org/gatk/). I have
   used versions in the 3.2 series and version 3.3-0. The specific jar file to use can be
   specified at the command line.

## Usage

If python2.7 is in your path, and if the script has the executable attribute, then you should
be able to execute `./regions-poorly-covered.py -h` to get a list of options. If python2.7 is not
in your executable path, then you may need to execute `/path/to/python2.7 regions-poorly-covered.py -h`

Required options:

* `-b <bamfiles.list>`

  A file containing the names of the bam files to analyze, one per line

* `-i <intervals>`

  A GATK-styled interval_list file or (0-based) bed file of intervals of interest

* `-o <output_basename>`

  The prefix to use for output files

* `-r <reference.fasta>`

  The reference genome sequence file against which to the bam files were aligned
  (GATK uses this to determine chromosome sequence lengths. The fasta file must
  have an associated ".dict" file per the GATK documentation)

Options with reasonable defaults:

* `-c <coverage>`

  The threshold for "poor coverage" (defaults to 20)

* `-p <percent>`

  Only output lines where percentage of samples with coverage < "-c" is above this value (defaults to 90)

* `-q <base_quality>`

  Minimum base quality to be considered for coverage calculations (--minBaseQuality in DoC; defaults to 20)

* `-Q <mapping_quality>`

  Minimum mapping quality to be considered for coverage calculations (--minMappingQuality in DoC; defaults to 20)

* `-m <memory>`

  Memory to allocate the the java process (used as java "-Xmx" option, defaults to 18g)

* `-t <threads>`

  Number of threads to use in DepthOfCoverage (defaults to 1, but should scale well to more threads...).
  
  **WARNING:** if the value for `threads` is set too high, then you may get a java exception about too many
  files being open. I use 4 threads for 200 samples on a cluster without difficulty.

* `-g <gatk.jar>`

  Path to the GATK '.jar' file (defaults to `GenomeAnalysisTK.jar` in the working directory)

* `-j <java_executable>`

  Defaults to the `java` executable in your current $PATH

Special options:

* `-s <previous_coverage_file>`

  Skip running DepthOfCoverage (the longest part of the whole task), by using an existing
  '.coverage' or '.coverage.gz' file. This could be useful if you want to re-filter with
  different values of `-c` or `-p` (if you want to change, `-q` or `-Q`, then DoC needs
  to be run with those different paramenters).

* `-f`

  Use a FIFO to allow the '.coverage' file to be compressed on-the-fly before being written
  to disk. This would be of benefit if IO is more limiting than CPU (which it often is,
  especially on network filesystems). If unspecified, the '.coverage' file will be
  written out at its full size

## An aside about the `-f` option

A FIFO is a `special file` that allows for a form of interprocess communication by having one
process write to the file and another read from it. Since DoC produces a very large (but
very compressable) output file that cannot be directly processed in a shell pipeline
(because it does not go to stdout or stderr), using a FIFO is the only way I could
think of to avoid a lot more (potentially network) disk IO. The `-f` option may not work if
FIFOs are not supported on a particular filesystem, in which case you should just try
running without it.

## output files

* Output files from DepthOfCoverage
* `*.DoC.c20.p90.txt` : the `*.coverage` file filtered when more than 90% of samples have coverage less than 20x
* `*.DoC.c20.p90.txt.bed` : a (0-based) bed file with the sites from the above file merged into continuous regions,
  also providing the average coverage across the region and the average % of samples with poor coverage across
  each region.

## TL;DR

Try running this (replacing things between `<angle-brackets>` with appropriate file names):

```
./regions-poorly-covered.py -b <bamfile.list> -r <reference.fasta> -i ACMG_IF_2014_inHGMD.0based.bed -o <basename>.acmg -c 20 -p 20 -q 20 -Q 20 -m 18g -t <threads> -g <GenomeAnalysisTK.jar> -f
```

Then, if that works, you could also try:
```
./regions-poorly-covered.py -b <bamfile.list> -r <reference.fasta> -i GeneTests_2015Feb23.0based.bed -o <basename>.genetests -c 20 -p 20 -q 20 -Q 20 -m 18g -t <threads> -g <GenomeAnalysisTK.jar> -f
```

Many of the options in both of these commands just specify the defaults, but they are spelled out here for clarity...

## License

The python script here is rather simple, I am making it available under a MIT-style license because I think
all open-source code should have some license, and this is one of the least-restrictive licenses out there.

GATK, of course, has its own license...
