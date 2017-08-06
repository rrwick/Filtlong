//Copyright 2017 Ryan Wick

//This file is part of LongQC

//LongQC is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//LongQC is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with LongQC.  If not, see <http://www.gnu.org/licenses/>.


#include <iostream>
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <limits>

#include "kseq.h"
#include "read.h"
#include "arguments.h"
#include "kmers.h"

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
    std::vector<Read> reads;
    if (!args.verbose)
        std::cerr << "Scoring long reads\n";
    int l;
    gzFile fp = gzopen(args.input_reads.c_str(), "r");
    kseq_t * seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if (l == -3)
            std::cerr << "Error reading " << args.input_reads << "\n";
        else {
            reads.emplace_back(seq->name.s, seq->seq.s, seq->qual.s, seq->seq.l, &kmers, &args);
            if (args.verbose)
                reads.back().print_verbose_read_info();
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    // If the user set thresholds using either --target_bases or --keep_percent, then we need to see which additional
    // reads should be labelled as failed.
    if (args.target_bases_set || args.keep_percent_set) {

        // Gather up pointers to the reads. If a read has been trimmed/split, it's these child reads which we use, not
        // the parent read.
        std::vector<Read *> read_pointers;
        for (size_t i = 0; i < reads.size(); ++i) {
            Read * read = &(reads[i]);
            if (read->child_reads.size() == 0) {
                read_pointers.push_back(read);
            }
            else {
                for (auto child : read->child_reads)
                    read_pointers.push_back(child);
            }
        }

        // See how many bases have already been passed.
        long long total_bases = 0;
        long long passed_bases = 0;
        for (auto read : read_pointers) {
            total_bases += read->m_length;
            if (read->m_passed)
                passed_bases += read->m_length;
        }
        std::cerr << "\n";
        std::cerr << read_pointers.size() << " reads (" << total_bases << " bp)\n";

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
        std::cerr << "Aiming to keep " << target_bases << " bp of reads\n";
        if (target_bases >= total_bases) {
            std::cerr << "Not enough reads to reach target\n";
        }
        else if (target_bases >= passed_bases) {
            std::cerr << "Reads already fall below target after filtering\n";
        }
        else {
            // Sort reads from best to worst.
            std::sort(read_pointers.begin(), read_pointers.end(),
                      [](const Read* a, const Read* b) {return a->m_final_score > b->m_final_score;});

            // Fail all reads after the threshold has been met.
            long long bases_so_far = 0;
            for (auto read : read_pointers) {
                if (bases_so_far > target_bases)
                    read->m_passed = false;
                else if (read->m_passed)
                    bases_so_far += read->m_length;
            }
            std::cerr << "Keeping " << bases_so_far << " bp\n";
        }
    }

    // Read through input reads again, this time outputting the keepers to stdout and ignoring the failures.
    // If in verbose mode, display a table as we go.
    if (!args.verbose)
        std::cerr << "\nOutputting passed long reads\n";
    fp = gzopen(args.input_reads.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if (l == -3)
            std::cerr << "Error reading " << args.input_reads << "\n";
        else {
//            seq->name.s;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    return 0;
}
