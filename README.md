<p align="center"><img src="misc/filtlong_logo.png" alt="Filtlong" width="450"></p>

Filtlong is a tool for filtering long reads by quality. It can take a large set of long reads and produce a smaller, higher quality subset. It uses both read length and sequence quality to choose which to output.



## Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Example commands (quick)](#example-commands-quick)
* [Example commands (detailed)](#example-commands-detailed)
* [Full usage](#full-usage)
* [Method](#method)
* [Read scoring](#read-scoring)
* [Trimming and splitting](#trimming-and-splitting)
* [Acknowledgements](#acknowledgements)
* [License](#license)



## Requirements

* Linux or macOS
* C++ compiler (C++11 support is required)
* zlib (usually included with Linux/macOS)



## Installation

Filtlong builds into a stand-alone executable:
```
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
bin/filtlong -h
```

You can then optionally copy Filtlong to a directory in your PATH:
```
cp bin/filtlong ~/.local/bin
```


## Example commands (quick)


### Without an external reference

```
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```

### With an external reference

```
filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 250 input.fastq.gz | gzip > output.fastq.gz
```


## Example commands (detailed)

These examples use a 1.3 Gbp read set that's part of a [barcoded 1D MinION run](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing). I assessed read quality by aligning the reads to a completed assembly using [minimap2](https://github.com/lh3/minimap2). Here is what the read length and identity distribution looks like before running Filtlong:
<p align="center"><img src="misc/example_commands_0_unfiltered.png" alt="unfiltered" width="400"></p>

### Without an external reference

When no external reference is provided, Filtlong judges read quality using the Phred quality scores in the FASTQ file.

```
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```

* `--min_length 1000`<br>
Discard any read which is shorter than 1 kbp.
* `--keep_percent 90`<br>
Throw out the worst 10% of reads. This is measured by bp, not by read count. So this option throws out the worst 10% of read bases.
* `--target_bases 500000000`<br>
If there are still more than 500 Mbp after throwing out reads under 1 kbp and the worst 10%, the remove the worst reads until only 500 Mbp remain. Useful for very large read sets.
* `input.fastq.gz`<br>
The input long reads to be filtered.
* `| gzip > output.fastq.gz`<br>
Filtlong outputs the filtered reads to stdout. Pipe to gzip to keep the file size down.

After running, Filtlong has cut the 1.3 Gbp down to a much better 500 Mbp subset:
<p align="center"><img src="misc/example_commands_1_without_reference.png" alt="without_reference" width="400"></p>


### With Illumina read reference

When an external reference is provided, Filtlong ignores the Phred quality scores and instead judges read quality using k-mer matches to the reference. This is a more accurate gauge of quality and enables a couple more options.

```
filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```

* `-1 illumina_1.fastq.gz -2 illumina_2.fastq.gz`<br>
These options allow Illumina reads to be used as an external reference. You can instead use `-a` to provide an assembly as a reference, but Illumina reads are preferable if they are available.

By using an external reference, Filtlong is better able to judge read quality. This shows in the resulting read stats, where almost all reads are now above 85% identity:
<p align="center"><img src="misc/example_commands_2_with_reference.png" alt="with_reference" width="400"></p>



### With trimming and splitting

When an external reference is provided, you can turn on read trimming and splitting to further increase read quality. See [Trimming and splitting](#trimming-and-splitting) for more information.

```
filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 250 input.fastq.gz | gzip > output.fastq.gz
```

* `--trim`<br>
Trim low-quality bases from the start and end. In this context 'low-quality' means bases which do not match a k-mer in the Illumina reads. This ensures the each read starts and ends with solid sequence.
* `--split 250`<br>
Split reads whenever 250 consequence bases fail to match a k-mer in the Illumina reads. This serves to remove very poor parts of reads while keeping the good parts. A lower value will split more aggressively and a higher value will be more conservative.

Trimming and splitting has further increased the output read identity:
<p align="center"><img src="misc/example_commands_3_trim_split.png" alt="trim_split" width="400"></p>


### Length priority

You can adjust the relative importance of Filtlong's read metrics. In this example, more weight is given to read length.

```
filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 1000 --length_weight 10 input.fastq.gz | gzip > output.fastq.gz
```

* `--length_weight 10`<br>
A length weight of 10 (instead of the default of 1) makes read length the most important factor when choosing the best reads.
* `--split 1000`<br>
By using a larger split value, Filtlong is less likely to split a read. This helps to keep the output reads on the long side.

The resulting 500 Mbp of reads are now mostly over 20 kbp, though some have mediocre identity:
<p align="center"><img src="misc/example_commands_4_length_priority.png" alt="length_priority" width="400"></p>


### Quality priority

You can adjust the relative importance of Filtlong's read metrics. In this example, more weight is given to read length.

```
filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 100 --mean_q_weight 10 input.fastq.gz | gzip > output.fastq.gz
```

* `--mean_q_weight 10`<br>
A mean quality weight of 10 (instead of the default of 1) makes mean read quality the most important factor when choosing the best reads.
* `--split 100`<br>
By using a smaller split value, Filtlong will split reads more often. This results in shorter reads but of higher quality.

The resulting 500 Mbp of reads are now very high identity, though most are under 25 kbp:
<p align="center"><img src="misc/example_commands_5_quality_priority.png" alt="quality_priority" width="400"></p>


## Full usage

```
usage: filtlong {OPTIONS} [input_reads]

Filtlong: a quality filtering tool for Nanopore and PacBio reads

positional arguments:
   input_reads                          Input long reads to be filtered

optional arguments:
   output thresholds:
      -t[int], --target_bases [int]        keep only the best reads up to this many total bases
      -p[float], --keep_percent [float]    keep only this fraction of the best reads
      --min_length [int]                   minimum length threshold
      --min_mean_q [float]                 minimum mean quality threshold
      --min_window_q [float]               minimum window quality threshold

   external references (if provided, read quality will be determined using these instead of from the
   Phred scores):
      -a[file], --assembly [file]          reference assembly in FASTA format
      -1[file], --illumina_1 [file]        reference Illumina reads in FASTQ format
      -2[file], --illumina_2 [file]        reference Illumina reads in FASTQ format

   score weights (control the relative contribution of each score to the final read score):
      --length_weight [float]              weight given to the length score (default: 1)
      --mean_q_weight [float]              weight given to the mean quality score (default: 1)
      --window_q_weight [float]            weight given to the window quality score (default: 1)

   read manipulation:
      --trim                               trim reads non-k-mer-matching bases from start/end of reads
      --split [split]                      split reads when this many bases lack a k-mer match

   other:
      --window_size [int]                  size of sliding window used when measuring window quality
                                           (default: 250)
      --verbose                            print a table with info for each read
      --version                            display the program version and quit

   -h, --help                           display this help menu

For more information, go to: https://github.com/rrwick/Filtlong
```


## Method

When run, Filtlong carries out the following steps:

1. If an external reference was provided, hash all of the reference's 16-mers.
  * If the reference is an assembly, then Filtlong simply hashes all 16-mers in the assembly.
  * If the reference is in Illumina reads, then the 16-mer has to be encountered a few times before it's hashed (to avoid hashing 16-mers that result from read errors).
2. Score each of the input reads.
  * Each read gets a final scores which is a function of its length, mean quality and window quality (see [Read scoring](#read-scoring) for more information).
  * If a read fails to meet any of the hard thresholds (`--min_length`, `--min_mean_q` or `--min_window_q`) then it is marked as 'fail' now.
  * If `--trim` or `--split` was used, then 'child' reads are made here (see [Trimming and splitting](#trimming-and-splitting) for more information). Each child read is scored using the same read scoring logic.
  * If `--verbose` was used, display detailed information about the read scoring.
3. Gather up all reads eligible for output. If neither `--trim` nor `--split` was used, this is simply the original set of reads. If `--trim` or `--split` was used, then the child reads replace the original reads.
4. If `--target_bases` and/or `--keep_percent` was used, sort the reads by quality and set an appropriate score threshold. Reads which fall below the threshold are marked as 'fail'.
  * If both `--target_bases` and `--keep_percent` are used, then the threshold is set to the more stringent of the two.
5. Output all reads which are not marked as 'fail' to stdout.
  * Reads are outputted in the same order as the input file (not in in quality-sorted order).


## Read scoring

Reads are scored based on three separate metrics: length, mean quality and window quality:

* __Length score__<br>
The length score is pretty simple: longer is better and [here is a graph of the score function](https://www.desmos.com/calculator/5m1lwd7fye). A read length of 5 kbp is considered mediocre and gets a score of 50. Shorter reads get a lower score and the score approaches 100 as the read length goes to infinity.

* __Mean quality__<br>
The mean quality is calculated in two different ways, depending on whether an external reference was used.
  * If an external reference was not used, the mean quality score is the average read identity, as indicated by the Phred quality scores. For example, consider a read where all the fastq quality characters are `+`. The qscores for each base are 10 which equates to a 90% chance of being correct. This read would then have a mean quality score of 90.
  * If an external reference was used, then Filtlong tallies up the 16-mers in the reference. Read bases are then scored as either 100 (contained in a 16-mer from the reference) or 0 (not contained in a 16-mer from the reference). The mean read score is a mean of all these base scores.

* __Mean quality score__<br>
Read mean qualities are converted to a z-score and scaled to the range 0-100 to make the mean quality score. This means that regardless of the actual mean quality distribution, the read with the worst mean quality will get a mean quality score of 0 and the read with the best mean quality will get a mean quality score of 100.

* __Window quality__<br>
The window quality is the mean quality for the lowest scoring window in the read. Each window's quality is calculated in the same way as the mean quality score (just for the window instead of for the whole read). The default window size is 250 but can be changed with `--window_size`.

* __Window quality score__<br>
The window quality score is the mean quality score scaled down by the window quality to mean quality ratio.

The read's final score is then a combination of the three scores:
* First, Filtlong takes the geometric weighted mean of the length score and the mean quality score. The weights are equal (both 1) by default, but you can adjust this with `--length_weight` and `--mean_q_weight` to make length or quality more or less important. By using a geometric mean instead of an arithmetic mean, the mean will be closer to the weaker of the two component scores.
* Then the score is scaled down by the ratio of the window quality score to the mean quality score.
* Expressed mathematically, the formula for the final read score is:
<p align="center"><img src="misc/score_equation.png" alt="score equation" width="500"></p>

It is the read's final score which is used to determine thresholds for the `--keep_percent` and `--target_bases` options.



## Trimming and splitting

If `--trim` or `--split` are used, then each read can result in one or more 'child' reads. It is these child reads which are considered in the final analysis. These options can only be used with an external reference (i.e. Filtlong will not trim or split based on Phred quality scores).


### Trim example

Consider the following example read. Bases which match a 16-mer to the reference are bold:
<p align="center"><img src="misc/trim_example_1.png" alt="Trim example 1" width="100%"></p>

If Filtlong is run with `--trim`, any non-matching bases at the start and end are removed to produce this child read:
<p align="center"><img src="misc/trim_example_2.png" alt="Trim example 2" width="100%"></p>

Only the child read (not the original read) is now eligible for outputting. Its length score will be a bit worse than the original read, but its mean quality score will be better.


### Split example

Using the same example read, if Filtlong was run with `--trim --split 20`, then in addition to trimming the ends, any run of non-matching bases 20 or longer will be removed. This results in two separate child reads:
<p align="center"><img src="misc/split_example_1.png" alt="Split example 1" width="100%"></p>

These child reads will have worse length scores than the original read but much better quality scores. Note that a split setting of 20 was used for this toy example but is very low for a real run of Filtlong – it would be quite aggressive and result in reads being split into very many pieces. A value more like 250 is more practical.

This is what you'd get if you ran the example read with `--trim --split 1`:
<p align="center"><img src="misc/split_example_2.png" alt="Split example 2" width="100%"></p>

Now _any_ run of non-matching bases is removed, regardless of length. The read is split into 3 child reads, each with perfect quality. Again, such a low setting would probably not be practical for a real read set.


### Real read example

Here's a real example of Filtlong's trimming and splitting. The output shown is what you'd see if you ran Filtlong with the `--verbose` option.

The input read is quite long and got a very good length score. Its mean quality score is decent but the window quality is zero, indicating that the read has windows with no bases matching the reference 16-mers:
```
bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template
            length = 117786     mean quality = 53.02      window quality =  0.00
```

The read's bad ranges are the coordinates non-matching start/end regions (because of `--trim`) or runs of 250 or more non-matching bases (because of `--split 250`). The child ranges are the inverse: bases that aren't in the bad ranges:
```
        bad ranges = 0-25, 71401-71745, 72393-72683, 72742-73049, 77279-77627, 78710-79055, 85575-85947, 86620-86877, 89397-89682, 91451-91782, 94415-94764, 96010-96306, 96604-96886, 98691-99176, 102349-102655, 102913-103286, 103488-103828, 106124-106397, 113277-113581, 117784-117786
      child ranges = 25-71401, 71745-72393, 72683-72742, 73049-77279, 77627-78710, 79055-85575, 85947-86620, 86877-89397, 89682-91451, 91782-94415, 94764-96010, 96306-96604, 96886-98691, 99176-102349, 102655-102913, 103286-103488, 103828-106124, 106397-113277, 113581-117784
```

A child read is made from each child range and assessed separately. The original read is no longer eligible for output – only the child reads are. Here are the first few of them:
```
bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_26-71401
            length = 71376      mean quality = 74.09      window quality =  6.80

bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_71746-72393
            length = 648        mean quality = 15.12      window quality =  6.40

bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_72684-72742
            length = 59         mean quality = 98.31      window quality = 98.31

bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_73050-77279
            length = 4230       mean quality = 25.11      window quality =  6.40

bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_77628-78710
            length = 1083       mean quality = 14.68      window quality =  6.40

bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_79056-85575
            length = 6520       mean quality = 23.27      window quality =  6.40

bf09f0e9-d27d-4a18-bced-f2536b62b3e5_Basecall_Alignment_template_85948-86620
            length = 673        mean quality = 18.28      window quality =  6.80
```

You can see the first child read is about 60% of the original read and despite being shorter, its higher quality will earn it a better final score. Many of the remaining pieces are very short and will therefore earn low scores. Depending on the output thresholds (like `--min_length` and `--keep_percent`) these may fail and will be discarded.



## Acknowledgements

I owe many thanks to [Kat Holt](https://holtlab.net/) and [Louise Judd](https://scholar.google.com.au/citations?user=eO22mYUAAAAJ) for keeping me well supplied with Nanopore reads.

Filtlong makes use of some nice open source libraries:
* [Klib](https://github.com/attractivechaos/klib/) for easy fastq parsing
* [args](https://github.com/Taywee/args) for command-line argument parsing.
* [C++ bloom filter library](https://github.com/ArashPartow/bloom) for some memory-saving in the k-mer counting

Thank you to the developers of these libraries!



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
