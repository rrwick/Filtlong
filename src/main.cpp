// Copyright 2017 Ryan Wick

// This file is part of LongQC

// LongQC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
// version.

// LongQC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.

// You should have received a copy of the GNU General Public License along with LongQC.  If not, see
// <http://www.gnu.org/licenses/>.


#include <iostream>
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <limits>
#include <unordered_map>
#include <utility>

#include "kseq.h"
#include "read.h"
#include "arguments.h"
#include "kmers.h"
#include "misc.h"

#define PROGRAM_VERSION "0.1.0"

KSEQ_INIT(gzFile, gzread)


int main(int argc, char **argv)
{
    Arguments args(argc, argv);
    if (args.parsing_result == BAD)
        return 1;
    else if (args.parsing_result == HELP)
        return 0;
    else if (args.parsing_result == VERSION) {
        std::cout << "LongQC v" << PROGRAM_VERSION << "\n";
        return 0;
    }

//    // TEMP - check argument parsing
//    std::cerr << "\n\n";
//    std::cerr << "input_reads: " << args.input_reads << std::endl;
//    if (args.target_bases_set) { std::cerr << "target_bases: " << args.target_bases << std::endl; }
//    else { std::cerr << "target_bases: not set" << std::endl; }
//    if (args.keep_percent_set) { std::cerr << "keep_percent: " << args.keep_percent << std::endl; }
//    else { std::cerr << "keep_percent: not set" << std::endl; }
//    if (args.assembly_set) { std::cerr << "assembly: " << args.assembly << std::endl; }
//    else { std::cerr << "assembly: not set" << std::endl; }
//    if (args.illumina_reads.size() > 0) {
//        std::cerr << "illumina_reads: " << std::endl;
//        for (auto i : args.illumina_reads)
//            std::cerr << "    " << i << std::endl;
//    }
//    else { std::cerr << "illumina_reads: not set" << std::endl; }
//    if (args.min_length_set) { std::cerr << "min_length: " << args.min_length << std::endl; }
//    else { std::cerr << "min_length: not set" << std::endl; }
//    if (args.min_mean_q_set) { std::cerr << "min_mean_q: " << args.min_mean_q << std::endl; }
//    else { std::cerr << "min_mean_q: not set" << std::endl; }
//    if (args.min_window_q_set) { std::cerr << "min_window_q: " << args.min_window_q << std::endl; }
//    else { std::cerr << "min_window_q: not set" << std::endl; }
//    std::cerr << "length_weight: " << args.length_weight << std::endl;
//    std::cerr << "mean_q_weight: " << args.mean_q_weight << std::endl;
//    std::cerr << "window_q_weight: " << args.window_q_weight << std::endl;
//    std::cerr << "trim: " << args.trim << std::endl;
//    if (args.split_set) { std::cerr << "split: " << args.split << std::endl; }
//    else { std::cerr << "split: not set" << std::endl; }
//    std::cerr << "window_size: " << args.window_size << std::endl;
//    std::cerr << "verbose: " << args.verbose << std::endl;
//    return 0;

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
    while ((l = kseq_read(seq)) >= 0) {
        if (l == -3)
            std::cerr << "Error reading " << args.input_reads << "\n";
        else {
            total_bases += seq->seq.l;
            std::string read_name = seq->name.s;
            Read * read = new Read(read_name, seq->seq.s, seq->qual.s, seq->seq.l, &kmers, &args);
            reads.push_back(read);
            if (args.verbose)
                read->print_verbose_read_info();

            if (read_dict.find(read->m_name) != read_dict.end()) {
                std::cerr << "Error: duplicate read name: " << read->m_name << "\n";
                return 0;
            }
            read_dict[read->m_name] = read;

            if (total_bases - last_progress >= 483611) {  // a big prime number so progress updates don't round off
                last_progress = total_bases;
                print_read_score_progress(reads.size(), total_bases);
            }
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    print_read_score_progress(reads.size(), total_bases);
    std::cerr << "\n";

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
                if (bases_so_far > target_bases)
                    read->m_passed = false;
                else if (read->m_passed)
                    bases_so_far += read->m_length;
            }
            std::cerr << "  keeping " << int_to_string(bases_so_far) << " bp\n";
        }
        std::cerr << "\n";
    }

    // Read through input reads again, this time outputting the keepers to stdout and ignoring the failures.
    // If in verbose mode, display a table as we go.
    if (!args.verbose)
        std::cerr << "Outputting passed long reads\n";
    fp = gzopen(args.input_reads.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if (l == -3)
            std::cerr << "Error reading " << args.input_reads << "\n";
        else {
            Read * read = read_dict[seq->name.s];

            if (read->m_child_reads.size() == 0) {
                if (read->m_passed) {
                    std::cout << "@" << seq->name.s;
                    if (seq->comment.l > 0)
                        std::cout << " " << seq->comment.s;
                    std::cout << "\n";
                    std::cout << seq->seq.s << "\n";
                    std::cout << "+\n";
                    std::cout << seq->qual.s << "\n";
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
                            std::cout << "@" << child_read->m_name;
                            if (seq->comment.l > 0)
                                std::cout << " " << seq->comment.s;
                            std::cout << "\n";

                            std::string seq_str = seq->seq.s;
                            std::string qual_str = seq->qual.s;
                            std::cout << seq_str.substr(start, length) << "\n";
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
