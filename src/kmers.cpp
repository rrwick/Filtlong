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


#include "kmers.h"

#include <iostream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


Kmers::Kmers() {

}


void Kmers::add_read_fastqs(std::vector<std::string> filenames) {
    std::cerr << "Hashing k-mers from Illumina reads\n";

    int sequence_count = 0;
    for (auto & filename : filenames) {
        std::cout << "  " << filename << "\n";
        sequence_count += add_reference(filename, true);
    }
    std::cout << "  " << m_first_time_kmers.size() << " total k-mers, ";
    std::cout << m_kmers.size() << " final k-mers\n";
}


void Kmers::add_assembly_fasta(std::string filename) {
    std::cerr << "Hashing k-mers from assembly\n";
    std::cout << "  " << filename << "\n";
    int sequence_count = add_reference(filename, false);
    std::cout << "  " << sequence_count << " contigs, ";
    std::cout << m_kmers.size() << " k-mers\n";
}


int Kmers::add_reference(std::string filename, bool require_two_kmer_copies) {
    int l;
    uint32_t forward_kmer, reverse_kmer;
    int sequence_count = 0;

    // We'll use a different k-mer adding function for assembly hashing and read hashing.
    void (Kmers::*add_kmer)(uint32_t);
    if (require_two_kmer_copies)
        add_kmer = &Kmers::add_kmer_require_two_copies;
    else
        add_kmer = &Kmers::add_kmer_require_one_copy;

    gzFile fp = gzopen(filename.c_str(), "r"); // STEP 2: open the file handler
    kseq_t * seq = kseq_init(fp);              // STEP 3: initialize seq
    while ((l = kseq_read(seq)) >= 0) {        // STEP 4: read sequence
        if (l == -3)
            std::cerr << "Error reading " << filename << "\n";
        else {
            ++sequence_count;
//            std::cout << "  name: " << seq->name.s;           // TEMP
//            if (seq->comment.l)                               // TEMP
//                std::cout << " " << seq->comment.s << "\n";   // TEMP
//            std::cout << "  length: " << seq->seq.l << "\n";  // TEMP

            if (seq->seq.l < 16)
                continue;

            char * sequence = seq->seq.s;

            // Build the starting k-mers from the first 16 bases.
            forward_kmer = starting_kmer_to_bits_forward(sequence);
            reverse_kmer = starting_kmer_to_bits_reverse(sequence);

            (this->*add_kmer)(forward_kmer);
            (this->*add_kmer)(reverse_kmer);

//            std::bitset<32> x(forward_kmer);                  // TEMP
//            std::bitset<32> y(reverse_kmer);                  // TEMP
//            std::cout << x << " " << y << "\n";               // TEMP

            for (size_t i = 16; i < seq->seq.l; ++i) {
                forward_kmer <<= 2;
                forward_kmer |= base_to_bits_forward(sequence[i]);

                reverse_kmer >>= 2;
                reverse_kmer |= base_to_bits_reverse(sequence[i]);

                (this->*add_kmer)(forward_kmer);
                (this->*add_kmer)(reverse_kmer);

//                std::bitset<32> x(forward_kmer);              // TEMP
//                std::bitset<32> y(reverse_kmer);              // TEMP
//                std::cout << x << " " << y << "\n";           // TEMP
            }
//            std::cout << "\n";                                  // TEMP
        }
    }
    kseq_destroy(seq);                         // STEP 5: destroy seq
    gzclose(fp);                               // STEP 6: close the file handler

    return sequence_count;
}


void Kmers::add_kmer_require_one_copy(uint32_t kmer) {
    m_kmers.insert(kmer);
}


void Kmers::add_kmer_require_two_copies(uint32_t kmer) {
    if (m_first_time_kmers.find(kmer) != m_first_time_kmers.end())  // if the k-mer has been seen before...
        m_kmers.insert(kmer);
    else
        m_first_time_kmers.insert(kmer);
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
