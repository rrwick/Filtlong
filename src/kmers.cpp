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


#include "kmers.h"

#include <iostream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "misc.h"

KSEQ_INIT(gzFile, gzread)


Kmers::Kmers() {
    bloom_parameters parameters;

    // TO DO: it might be worth experimenting with these values to see how it affects time and memory usage.
    parameters.projected_element_count = 100000000;
    parameters.false_positive_probability = 0.0001; // 1 in 10000
    parameters.random_seed = 0xA5A5A5A5;

    parameters.compute_optimal_parameters();

    //Instantiate Bloom Filter
    bloom = new bloom_filter(parameters);

    required_kmer_copies = 4;
}


Kmers::~Kmers() {
    delete bloom;
}


void Kmers::add_read_fastqs(std::vector<std::string> filenames) {
    std::cerr << "Hashing 16-mers from Illumina reads\n";

    int sequence_count = 0;
    for (auto & filename : filenames)
        sequence_count += add_reference(filename, true);
    std::cerr << "  " << int_to_string(sequence_count) << " reads, "
              << int_to_string(m_kmers.size()) << " 16-mers\n\n";
}


void Kmers::add_assembly_fasta(std::string filename) {
    std::cerr << "Hashing 16-mers from assembly\n";
    std::cerr << "  " << filename << "\n";
    int sequence_count = add_reference(filename, false);
    std::string noun;
    if (sequence_count == 1)
        noun = "contig";
    else
        noun = "contigs";
    std::cerr << "  " << int_to_string(sequence_count) << " " << noun << ", "
              << int_to_string(m_kmers.size()) << " 16-mers\n\n";
}


int Kmers::add_reference(std::string filename, bool require_two_kmer_copies) {
    int l;
    uint32_t forward_kmer, reverse_kmer;
    int sequence_count = 0;

    // We'll use a different k-mer adding function for assembly hashing and read hashing.
    void (Kmers::*add_kmer)(uint32_t);
    if (require_two_kmer_copies)
        add_kmer = &Kmers::add_kmer_require_multiple_copies;
    else
        add_kmer = &Kmers::add_kmer_require_one_copy;

    long long base_count = 0;
    long long last_progress = 0;

    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t * seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if (l == -3)
            std::cerr << "Error reading " << filename << "\n";
        else {
            ++sequence_count;

            // Can't get a 16-mer from a sequence shorter than 16 bp.
            if (seq->seq.l < 16)
                continue;

            base_count += seq->seq.l;
            char * sequence = seq->seq.s;

            // Build the starting k-mers from the first 16 bases.
            forward_kmer = starting_kmer_to_bits_forward(sequence);
            reverse_kmer = starting_kmer_to_bits_reverse(sequence);

            (this->*add_kmer)(forward_kmer);
            (this->*add_kmer)(reverse_kmer);

            for (size_t i = 16; i < seq->seq.l; ++i) {
                forward_kmer <<= 2;
                forward_kmer |= base_to_bits_forward(sequence[i]);

                reverse_kmer >>= 2;
                reverse_kmer |= base_to_bits_reverse(sequence[i]);

                (this->*add_kmer)(forward_kmer);
                (this->*add_kmer)(reverse_kmer);
            }

            if (base_count - last_progress >= 483611) {  // a big prime number so progress updates don't round off
                last_progress = base_count;
                print_hash_progress(filename, base_count);
            }
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    print_hash_progress(filename, base_count);
    std::cerr << "\n";
    return sequence_count;
}


void Kmers::add_kmer_require_one_copy(uint32_t kmer) {
    m_kmers.insert(kmer);
}


void Kmers::add_kmer_require_multiple_copies(uint32_t kmer) {
    // If the kmer is already in the final set, then we can skip the rest of this function.
    if (m_kmers.find(kmer) != m_kmers.end())
        return;

    // Check the bloom filter. If it's not in there, this is definitely the first time it's been seen.
    if (!bloom->contains(kmer))
        bloom->insert(kmer);

    // If it's in the bloom filter, then it's probably been seen once before (though maybe not, based on the false
    // positive rate of the bloom filter. Next we check the k-mer counts. If it's not in there, we say it's the second
    // time the kmer's been seen.
    else if (m_kmer_counts.find(kmer) == m_kmer_counts.end())
        m_kmer_counts[kmer] = 2;

    // If the k-mer is in the counts, then we increment its count. If the count is high enough, we add it to the k-mer
    // set (and remove it from the counts to save some memory).
    else {
        int times_seen = ++m_kmer_counts[kmer];
        if (times_seen >= required_kmer_copies) {
            m_kmers.insert(kmer);
            m_kmer_counts.erase(kmer);
        }
    }
}



bool Kmers::is_kmer_present(uint32_t kmer) {
    return m_kmers.find(kmer) != m_kmers.end();
}



uint32_t Kmers::base_to_bits_forward(char base) {
    switch (base) {
        case 'A':
            return 0;  // 00000000000000000000000000000000
        case 'C':
            return 1;  // 00000000000000000000000000000001
        case 'G':
            return 2;  // 00000000000000000000000000000010
        case 'T':
            return 3;  // 00000000000000000000000000000011
        case 'a':
            return 0;
        case 'c':
            return 1;
        case 'g':
            return 2;
        case 't':
            return 3;
    }
    return 0;
}


uint32_t Kmers::base_to_bits_reverse(char base) {
    switch (base) {
        case 'T':
            return 0;           // 00000000000000000000000000000000
        case 'G':
            return 1073741824;  // 01000000000000000000000000000000
        case 'C':
            return 2147483648;  // 10000000000000000000000000000000
        case 'A':
            return 3221225472;  // 11000000000000000000000000000000
        case 't':
            return 0;
        case 'g':
            return 1073741824;
        case 'c':
            return 2147483648;
        case 'a':
            return 3221225472;
    }
    return 0;
}


uint32_t Kmers::starting_kmer_to_bits_forward(char * sequence) {
    uint32_t kmer = 0;
    for (int i = 0; i < 16; ++i) {
        kmer <<= 2;
        kmer |= base_to_bits_forward(sequence[i]);
    }
    return kmer;
}


uint32_t Kmers::starting_kmer_to_bits_reverse(char * sequence) {
    uint32_t kmer = 0;
    for (int i = 0; i < 16; ++i) {
        kmer >>= 2;
        kmer |= base_to_bits_reverse(sequence[i]);
    }
    return kmer;
}
