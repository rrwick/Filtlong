#!/usr/bin/env bash

long_reads=$1
illumina_1=$2
illumina_2=$3

# Get the directory of this script.
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get paths to other necessary files.
hist="$script_dir"/histogram.py
filtlong="$script_dir"/../bin/filtlong

# Prepare a tab-delimited file of read info.
if [ -n "$illumina_1" ] && [ -n "$illumina_2" ]; then
    illumina_args="-1 "$illumina_1" -2 "$illumina_2
else
    illumina_args=""
fi
$filtlong $illumina_args --min_mean_q 1 --verbose $long_reads 2>&1 >/dev/null | grep -P "length = \d+"  | sed 's/\s\+length = //' | sed 's/\s\+\w\+ quality = /\t/g' > temp_filtlong_read_info

printf "\n"
echo "READ SET SUMMARY"
echo "----------------"
printf "number of reads: %'d\n" $(cat temp_filtlong_read_info | wc -l)
total_bases=$(cut -f1 temp_filtlong_read_info | paste -sd+ - | bc)
printf "number of bases: %'d\n" $total_bases
target=$(($total_bases / 2))
n50=$(zless $long_reads | paste - - - - | cut -f2 | awk '{print length($0);}' | sort -nr | awk -v tar="$target" '{sum += $0; if (sum > tar) {print($0); exit;}}')
printf "N50 read length: %'d\n" $n50

printf "\n\n"
echo "READ LENGTHS"
echo "------------"
cat temp_filtlong_read_info | awk '{print $1, $1}' | $hist -a -m 0 -b 40

printf "\n\n"
echo "MEAN QUALITIES"
echo "--------------"
cat temp_filtlong_read_info | awk '{print $1, $2}' | $hist -a -b 40

printf "\n\n"
echo "WINDOW QUALITIES"
echo "----------------"
cat temp_filtlong_read_info | awk '{print $1, $3}' | $hist -a -b 40

printf "\n"
rm temp_filtlong_read_info
