"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Filtlong

This module contains some tests for Filtlong. To run them, execute `python3 -m unittest` from the
root Filtlong directory.

This file is part of Filtlong. Filtlong is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Filtlong is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Filtlong. If
not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import os
import subprocess


def load_fastq(filename):
    reads = []
    with open(filename, 'rb') as fastq:
        for line in fastq:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            if not stripped_line.startswith(b'@'):
                continue
            name = stripped_line[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads.append((name, sequence, qualities))
    return reads


class TestSplit(unittest.TestCase):

    def run_command(self, command):
        binary_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'bin', 'filtlong')
        input_path = os.path.join(os.path.dirname(__file__), 'test_split.fastq')
        assembly_reference = os.path.join(os.path.dirname(__file__), 'test_reference.fasta')
        illumina_reference_1 = os.path.join(os.path.dirname(__file__), 'test_reference_1.fastq.gz')
        illumina_reference_2 = os.path.join(os.path.dirname(__file__), 'test_reference_2.fastq.gz')

        command = command.replace('filtlong', binary_path)
        command = command.replace('INPUT', input_path)
        command = command.replace('ASSEMBLY', assembly_reference)
        command = command.replace('ILLUMINA_1', illumina_reference_1)
        command = command.replace('ILLUMINA_2', illumina_reference_2)

        output_name = 'TEMP_' + str(os.getpid())
        command = command.replace('OUTPUT', output_name)
        try:
            self.output_file = [x for x in command.split() if output_name in x][0]
        except IndexError:
            self.output_file = ''
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, err = p.communicate()
        return err.decode()

    def tearDown(self):
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

    def check_one_read(self, read, length, start_seq, end_seq):
        self.assertEqual(len(read[1]), length)
        self.assertEqual(len(read[2]), length)
        self.assertTrue(read[1].startswith(start_seq))
        self.assertTrue(read[1].endswith(end_seq))

    def test_split_1(self):
        """
        When splitting is off, the reads aren't split.
        """
        console_out = self.run_command('filtlong --min_length 1 -a ASSEMBLY INPUT > OUTPUT.fastq')
        self.assertTrue('4 reads (11,600 bp)' in console_out)
        self.assertTrue('after splitting' not in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 4)

    def test_split_2(self):
        """
        When the split threshold is over 200, no reads are split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 250 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 4 reads (11,600 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 4)
        for i in range(0, 4):
            self.check_one_read(split_reads[i], 2900, b'TCATTACG', b'AAAAGGAC')

    def test_split_3(self):
        """
        When the split threshold is over 200, no reads are split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 201 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 4 reads (11,600 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 4)
        for i in range(0, 4):
            self.check_one_read(split_reads[i], 2900, b'TCATTACG', b'AAAAGGAC')

    def test_split_4(self):
        """
        When the split threshold is exactly 200, only the last read is split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 200 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 5 reads (11,400 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 5)
        for i in range(0, 3):
            self.check_one_read(split_reads[i], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[3], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[4], 1700, b'CCATGACA', b'AAAAGGAC')

    def test_split_5(self):
        """
        When the split threshold is between 150 and 200, only the last read is split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 175 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 5 reads (11,400 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 5)
        for i in range(0, 3):
            self.check_one_read(split_reads[i], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[3], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[4], 1700, b'CCATGACA', b'AAAAGGAC')

    def test_split_6(self):
        """
        When the split threshold is between 50 and 100, tbe last two reads are split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 75 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 6 reads (11,300 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 6)
        for i in range(0, 2):
            self.check_one_read(split_reads[i], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[2], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[3], 1800, b'GAAACGTG', b'AAAAGGAC')
        self.check_one_read(split_reads[4], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[5], 1700, b'CCATGACA', b'AAAAGGAC')

    def test_split_7(self):
        """
        When the split threshold is between 50 and 100, tbe last two reads are split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 51 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 6 reads (11,300 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 6)
        for i in range(0, 2):
            self.check_one_read(split_reads[i], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[2], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[3], 1800, b'GAAACGTG', b'AAAAGGAC')
        self.check_one_read(split_reads[4], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[5], 1700, b'CCATGACA', b'AAAAGGAC')

    def test_split_8(self):
        """
        When the split threshold is exactly 50, tbe last three reads are split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 50 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 7 reads (11,250 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 7)
        self.check_one_read(split_reads[0], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[1], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[2], 1850, b'TGATGAAT', b'AAAAGGAC')
        self.check_one_read(split_reads[3], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[4], 1800, b'GAAACGTG', b'AAAAGGAC')
        self.check_one_read(split_reads[5], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[6], 1700, b'CCATGACA', b'AAAAGGAC')

    def test_split_9(self):
        """
        When the split threshold is less than 50, tbe last three reads are split.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 25 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 7 reads (11,250 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 7)
        self.check_one_read(split_reads[0], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[1], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[2], 1850, b'TGATGAAT', b'AAAAGGAC')
        self.check_one_read(split_reads[3], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[4], 1800, b'GAAACGTG', b'AAAAGGAC')
        self.check_one_read(split_reads[5], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[6], 1700, b'CCATGACA', b'AAAAGGAC')
        
    def test_split_10(self):
        """
        Same as the previous test, but using reads as a reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --split 25 INPUT > OUTPUT.fastq')
        self.assertTrue('after splitting: 7 reads (11,250 bp)' in console_out)
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 7)
        self.check_one_read(split_reads[0], 2900, b'TCATTACG', b'AAAAGGAC')
        self.check_one_read(split_reads[1], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[2], 1850, b'TGATGAAT', b'AAAAGGAC')
        self.check_one_read(split_reads[3], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[4], 1800, b'GAAACGTG', b'AAAAGGAC')
        self.check_one_read(split_reads[5], 1000, b'TCATTACG', b'ATTGGTTG')
        self.check_one_read(split_reads[6], 1700, b'CCATGACA', b'AAAAGGAC')

    def test_split_names(self):
        """
        Makes sure the split reads are correctly named.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --split 25 INPUT > OUTPUT.fastq')
        split_reads = load_fastq(self.output_file)
        self.assertEqual(len(split_reads), 7)
        self.assertEqual(split_reads[0][0], b'test_split_1')
        self.assertEqual(split_reads[1][0], b'test_split_2_1-1000')
        self.assertEqual(split_reads[2][0], b'test_split_2_1051-2900')
        self.assertEqual(split_reads[3][0], b'test_split_3_1-1000')
        self.assertEqual(split_reads[4][0], b'test_split_3_1101-2900')
        self.assertEqual(split_reads[5][0], b'test_split_4_1-1000')
        self.assertEqual(split_reads[6][0], b'test_split_4_1201-2900')
