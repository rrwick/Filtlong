// Copyright 2017 Ryan Wick

// This file is part of Filtlong

// Filtlong is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
// version.

// Filtlong is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.

// You should have received a copy of the GNU General Public License along with Filtlong.  If not, see
// <http://www.gnu.org/licenses/>.


#include <iostream>
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <limits>
#include <unordered_map>
#include <utility>
#include <math.h>

#include "kseq.h"
#include "read.h"
#include "arguments.h"
#include "kmers.h"
#include "misc.h"

#define PROGRAM_VERSION "0.2.1"

KSEQ_INIT(gzFile, gzread)


int main(int argc, char **argv)
{
    Arguments args(argc, argv);
    if (args.parsing_result == BAD)
        return 1;
    else if (args.parsing_result == HELP)
        return 0;
    else if (args.parsing_result == VERSION) {
        std::cout << "Filtlong v" << PROGRAM_VERSION << "\n";
        return 0;
    }

    std::cerr << "\n";

    // Read through references and save 16-mers. For assembly references, this will save all 16-mers in the assembly.
    // For Illumina read references, the k-mer needs to appear a few times before it's added to the set.
    Kmers kmers;
    if (args.assembly_set || args.illumina_reads.size() > 0) {
        if (args.assembly_set)
            kmers.add_assembly_fasta(args.assembly);
        if (args.illumina_reads.size() > 0)
            kmers.add_read_fastqs(args.illumina_reads);
    }

    // Read through input long reads once, storing them as Read objects and calculating their scores.
    // While we go, make sure there are no duplicate read names. Quit with an error if so.
    long long total_bases = 0;
    long long last_progress = 0;
    std::vector<Read*> reads;
    std::unordered_map<std::string, Read*> read_dict;
    if (!args.verbose)
        std::cerr << "Scoring long reads\n";
    int l;
    gzFile fp = gzopen(args.input_reads.c_str(), "r");
    kseq_t * seq = kseq_init(fp);

    bool any_fasta = false;
    bool any_fastq = false;

    while (true) {
        l = kseq_read(seq);
        if (l == -1)  // end of file
            break;
        if (l == -2) {
            std::cerr << "Error: incorrect FASTQ format for read " << seq->name.s << "\n";
            return 1;
        }
        if (l == -3) {
            std::cerr << "Error reading " << args.input_reads << "\n";
            return 1;
        }
        else {
            total_bases += seq->seq.l;
            std::string read_name = seq->name.s;

            bool fasta_format = (seq->qual.l == 0 && seq->seq.l > 0);
            bool fastq_format = (seq->qual.l > 0 && seq->seq.l > 0 && seq->qual.l == seq->seq.l);

            any_fasta = (any_fasta || fasta_format);
            any_fastq = (any_fastq || fastq_format);
            if (any_fasta && any_fastq) {
                std::cerr << "\n\n" << "Error: could not parse input reads" << "\n";
                std::cerr << "  problem occurred at read " << read_name << "\n";
                return 1;
            }

            if (fasta_format && kmers.empty()) {
                std::cerr << "\n\n" << "Error: FASTA input not supported without an external reference" << "\n";
                return 1;
            }

            Read * read = new Read(read_name, seq->seq.s, seq->qual.s, int(seq->seq.l), &kmers, &args);
            reads.push_back(read);
            if (args.verbose)
                read->print_verbose_read_info();

            if (read_dict.find(read->m_name) != read_dict.end()) {
                std::cerr << "Error: duplicate read name: " << read->m_name << "\n";
                return 1;
            }
            read_dict[read->m_name] = read;

            if (total_bases - last_progress >= 483611) {  // a big prime number so progress updates don't round off
                last_progress = total_bases;
                if (!args.verbose)
                    print_read_score_progress(reads.size(), total_bases);
            }
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    if (!args.verbose)
        print_read_score_progress(reads.size(), total_bases);
    std::cerr << "\n";

    // Determine the output format.
    bool fasta_output = any_fasta;
    bool fastq_output = any_fastq;

    // Gather up reads to output. If a read has been trimmed/split, it's these child reads which we use, not the
    // parent read.
    std::vector<Read*> reads2;
    for (auto read : reads) {
        if (read->m_child_reads.size() == 0) {
            reads2.push_back(read);
        }
        else {
            for (auto child : read->m_child_reads)
                reads2.push_back(child);
        }
    }
    size_t longest_read_name = 0;
    for (auto read : reads2) {
        if (read->m_name.size() > longest_read_name)
            longest_read_name = read->m_name.size();
    }

    // If --trim or --split was used, display some summary info here.
    if (args.trim || args.split_set) {
        long long total_after_trim_split = 0;
        for (auto read : reads2)
            total_after_trim_split += read->m_length;
        if (args.trim && args.split_set)
            std::cerr << "  after trimming and splitting: ";
        else if (args.trim)
            std::cerr << "  after trimming: ";
        else
            std::cerr << "  after splitting: ";
        std::cerr << int_to_string(reads2.size()) << " reads (" << int_to_string(total_after_trim_split) << " bp)\n";
    }
    std::cerr << "\n";

    // Go through the mean quality scores and find the min, max, mean and standard deviation.
    double min_quality = 100.0;
    double max_quality = 0.0;
    double quality_sum = 0.0;
    for (auto read : reads2) {
        quality_sum += read->m_mean_quality;
        if (read->m_mean_quality > max_quality)
            max_quality = read->m_mean_quality;
        if (read->m_mean_quality < min_quality)
            min_quality = read->m_mean_quality;
    }
    double mean_quality = quality_sum / reads2.size();
    double stdev_sum = 0.0;
    for (auto read : reads2) {
        double mean_diff = read->m_mean_quality - mean_quality;
        stdev_sum += mean_diff * mean_diff;
    }
    double stdev_quality = sqrt(stdev_sum / reads2.size());
    double min_z_score, max_z_score;
    if (stdev_quality > 0.0) {
        min_z_score = (min_quality - mean_quality) / stdev_quality;
        max_z_score = (max_quality - mean_quality) / stdev_quality;
    }
    else {
        min_z_score = 1.0;
        max_z_score = 1.0;
    }
    double max_min_z_diff = max_z_score - min_z_score;

    // Now normalise each read's quality scores.
    if (args.verbose)
        std::cerr << "\n\n" << "Read name" << "\t" << "Length score" << "\t" << "Mean quality score" << "\t"
                  << "Window quality score" << "\t" << "Final score" << "\n";
    for (auto read : reads2) {
        double window_ratio = read->m_window_quality / read->m_mean_quality;
        if (window_ratio > 1.0)
            window_ratio = 1.0;
        double quality_z_score = (read->m_mean_quality - mean_quality) / stdev_quality;
        read->m_mean_quality = 100.0 * (quality_z_score - min_z_score) / max_min_z_diff;
        read->m_window_quality = read->m_mean_quality * window_ratio;
        read->set_final_score(args.length_weight, args.mean_q_weight, args.window_q_weight);
        if (args.verbose)
            read->print_scores(longest_read_name);
    }
    if (args.verbose)
        std::cerr << "\n";

    // If the user set thresholds using either --target_bases or --keep_percent, then we need to see which additional
    // reads should be labelled as failed.
    if (args.target_bases_set || args.keep_percent_set) {
        std::cerr << "Filtering long reads\n";

        // See how many bases have already been passed.
        long long passed_bases = 0;
        for (auto read : reads2) {
            if (read->m_passed)
                passed_bases += read->m_length;
        }

        // Determine how many bases we should keep.
        long long target_bases;
        if (args.target_bases_set)
            target_bases = args.target_bases;
        else
            target_bases = std::numeric_limits<long long>::max();
        if (args.keep_percent_set) {
            long long keep_target = (long long)((args.keep_percent / 100.0) * total_bases);
            target_bases = std::min(target_bases, keep_target);
        }
        std::cerr << "  target: " << int_to_string(target_bases) << " bp\n";
        if (target_bases >= total_bases) {
            std::cerr << "  not enough reads to reach target\n";
        }
        else if (target_bases >= passed_bases) {
            std::cerr << "  reads already fall below target after filtering\n";
        }
        else {
            // Sort reads from best to worst.
            std::sort(reads2.begin(), reads2.end(),
                      [](const Read* a, const Read* b) {return a->m_final_score > b->m_final_score;});

            // Fail all reads after the threshold has been met.
            long long bases_so_far = 0;
            for (auto read : reads2) {
                if (read->m_passed && bases_so_far < target_bases)
                    bases_so_far += read->m_length;
                else
                    read->m_passed = false;
            }
            std::cerr << "  keeping " << int_to_string(bases_so_far) << " bp\n";
        }
        std::cerr << "\n";
    }

    // Read through input reads again, this time outputting the keepers to stdout and ignoring the failures.
    std::cerr << "Outputting passed long reads\n";
    fp = gzopen(args.input_reads.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        Read * read = read_dict[seq->name.s];

        if (read->m_child_reads.size() == 0) {
            if (read->m_passed) {
                std::cout << (fasta_output ? ">" : "@");
                std::cout << seq->name.s;
                if (seq->comment.l > 0)
                    std::cout << " " << seq->comment.s;
                std::cout << "\n";
                std::cout << seq->seq.s << "\n";
                if (fastq_output) {
                    std::cout << "+\n";
                    std::cout << seq->qual.s << "\n";
                }
            }
        }
        else {
            for (size_t i = 0; i < read->m_child_reads.size(); ++i) {
                Read * child_read = read->m_child_reads[i];
                if (child_read->m_passed) {
                    std::pair<int,int> child_read_range = read->m_child_read_ranges[i];
                    int start = child_read_range.first;
                    int end = child_read_range.second;
                    int length = end - start;
                    if (length > 0) {
                        std::cout << (fasta_output ? ">" : "@");
                        std::cout << child_read->m_name;
                        if (seq->comment.l > 0)
                            std::cout << " " << seq->comment.s;
                        std::cout << "\n";

                        std::string seq_str = seq->seq.s;
                        std::cout << seq_str.substr(start, length) << "\n";

                        if (fastq_output) {
                            std::string qual_str = seq->qual.s;
                            std::cout << "+\n";
                            std::cout << qual_str.substr(start, length) << "\n";
                        }
                    }
                }
            }
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    // Clean up.
    for (auto read : reads)
        delete read;

    std::cerr << "\n";
    return 0;
}
