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

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


Kmers::Kmers() {

}


void Kmers::add_read_fastq(std::string filename) {
    add_reference(filename, true);
}


void Kmers::add_assembly_fasta(std::string filename) {
    add_reference(filename, false);
}


bool Kmers::is_kmer_present(uint32_t kmer) {
    return m_kmers.find(kmer) != m_kmers.end();
}


void Kmers::add_reference(std::string filename, bool require_two_kmer_copies) {

    gzFile fp;
    kseq_t *seq;
    int l;

    fp = gzopen(filename.c_str(), "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        printf("name: %s\n", seq->name.s);
        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        printf("seq: %s\n", seq->seq.s);
        if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
    }
    printf("return value: %d\n", l);
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler

}
