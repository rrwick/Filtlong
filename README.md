# Filtlong

Filtlong is a tool for distilling a large set of long reads to a smaller, higher quality subset. It scores reads based on their length and quality to choose which to output. If an external reference is available, it can use that to more accurately assess read quality.



## Requirements

* Linux or macOS
* C++ compiler (C++11 support is required)
* zlib (usually included with Linux/macOS)



## Installation

Filtlong should be simple to build:
```
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
bin/longqc -h
```

You can then optionally copy Filtlong to a directory in your path:
```
cp bin/longqc /usr/local/bin
```


## Example command (without an external reference)

When no external reference is provided, Filtlong judges read quality using the Phred quality scores in the FASTQ file.

```
longqc --min_length 1000 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```

Explanation:
* `--min_length 1000`<br>
Discard any read which is shorter than 1 kbp. Since in this example we have Illumina reads, there's no need for long reads that aren't actually long.
* `--keep_percent 90`<br>
Throw out the worst 10% of reads.
* `--target_bases 500000000`<br>
If there are still more than 500 Mbp after throwing out reads under 1 kbp and the worst 10%, the remove the worst reads until only 500 Mbp remain. Useful for very large read sets.
* `input.fastq.gz`<br>
The input long reads to be filtered.
* `| gzip > output.fastq.gz`
Filtlong outputs the filtered reads to stdout. I just pipe to gzip to keep the file size down.



## Example command (with Illumina read reference)

When an external reference is provided, Filtlong ignores the Phred quality scores and instead judges read quality using k-mer matches to the reference. This is a more accurate gauge of quality and enables a couple more options.

```
longqc -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 250 input.fastq.gz | gzip > output.fastq.gz
```

Explanation:
* `-1 illumina_1.fastq.gz -2 illumina_2.fastq.gz`<br>
Illumina reads as an external reference. You can instead use `-a` to provide an assembly as a reference (but Illumina reads are preferable if they are available).
* `--trim`<br>
Trim low-quality bases from the start and end. In this context 'low-quality' means bases which do not match a k-mer in the Illumina reads. This ensures the each read starts and ends with solid sequence.
* `--split 250`<br>
Split reads whenever 250 consequence bases fail to match a k-mer in the Illumina reads. This serves to remove very poor parts of reads while keeping the good parts. A lower value will split more aggressively and a higher value will be more conservative.



## Full usage

```
usage: longqc {OPTIONS} [input_reads]

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


## Read scoring

Reads are scored based on three separate metrics: length, mean quality and window quality:

* __Length score__<br>
The length score is pretty simple: longer is better and [here is a graph of the score function](https://www.desmos.com/calculator/5m1lwd7fye). A read length of 5 kbp is considered mediocre and gets a score of 50. Shorter reads get a lower score and the score approaches 100 as the read length goes to infinity.

* __Mean quality score__<br>
The mean quality score is calculated in two different ways, depending on whether an external reference was used.
  * If an external reference was not used, the mean quality score is the average read identity, as indicated by the Phred quality scores. For example, consider a read where all the fastq quality characters are `+`. The qscores for each base are 10 which equates to a 90% chance of being correct. This read would then have a mean quality score of 90.
  * If an external reference was used, then Filtlong tallies up the 16-mers in the reference. Read bases are then scored as either 100 (contained in a 16-mer from the reference) or 0 (not contained in a 16-mer from the reference). The mean read score is a mean of all these base scores.

* __Window quality score__<br>
The window quality score is the mean quality score for the lowest scoring window in the read. Each window's quality is calcuated in the same way as the mean quality score (just for the window instead of for the whole read). The default window size is 250 but can be changed with `--window_size`.

The read's final score is then a combination of those three scores:
* First, Filtlong takes the geometric weighted mean of the length score and the mean quality score. The weights are equal (both 1) by default, but you can adjust this with `--length_weight` and `--mean_q_weight` to make length or quality more or less important. By using a geometric mean instead of an arithmetic mean, the mean will be closer to the weaker of the two component scores.
* Then the score is scaled down by the ratio of the window quality score to the mean quality score.
* Expressed mathematically, the formula for the final read score is:
<p align="center"><img src="misc/score_equation.png" alt="score equation" width="400"></p>

It is the read's final score which is used to determine thresholds for the `--keep_percent` and `--target_bases` options.



## Trimming and splitting

If `--trim` or `--split` are used, then each read can result in one or more 'child' reads. It is these child reads which are considered in the final analysis.

Trim example: TO DO

Split example: TO DO



# Acknowledgements

I owe many thanks to [Kat Holt](https://holtlab.net/) and [Louise Judd](https://scholar.google.com.au/citations?user=eO22mYUAAAAJ) for keeping me well supplied with Nanopore reads.

Filtlong makes use of some nice open source libraries:
* [Klib](https://github.com/attractivechaos/klib/) for easy fastq parsing
* [args](https://github.com/Taywee/args) for command-line argument parsing.
* [C++ bloom filter library](https://github.com/ArashPartow/bloom) for some memory-saving in the k-mer counting

Thank you to the developers of these libraries!



# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
