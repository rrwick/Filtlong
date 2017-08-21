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


class TestSort(unittest.TestCase):

    def run_command(self, command):
        binary_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'bin', 'filtlong')
        input_path = os.path.join(os.path.dirname(__file__), 'test_sort.fastq')
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

    def test_sort_1(self):
        """
        With a high input threshold, all three reads should be outputted.
        """
        console_out = self.run_command('filtlong --target_bases 100000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 100,000 bp' in console_out)
        self.assertTrue('not enough reads to reach target' in console_out)

    def test_sort_2(self):
        """
        With an input threshold of just over 10000, all three reads should be outputted.
        """
        console_out = self.run_command('filtlong --target_bases 10001 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 10,001 bp' in console_out)
        self.assertTrue('keeping 15,000 bp' in console_out)

    def test_sort_3(self):
        """
        With an input threshold of exactly 10000, only two reads should be outputted. Read 1 is the worst (by Phred
        score), so just reads 2 and 3.
        """
        console_out = self.run_command('filtlong --target_bases 10000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 10,000 bp' in console_out)
        self.assertTrue('keeping 10,000 bp' in console_out)

    def test_sort_4(self):
        """
        With an input threshold of just over 5000, two reads should be outputted.
        """
        console_out = self.run_command('filtlong --target_bases 5001 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 5,001 bp' in console_out)
        self.assertTrue('keeping 10,000 bp' in console_out)

    def test_sort_5(self):
        """
        With an input threshold of exactly 5000, only one read should be outputted: read 2 (best by Phred score).
        """
        console_out = self.run_command('filtlong --target_bases 5000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_2'])
        self.assertTrue('target: 5,000 bp' in console_out)
        self.assertTrue('keeping 5,000 bp' in console_out)

    def test_sort_6(self):
        """
        With an input threshold of 1, only one read should be outputted: read 2 (best by Phred score).
        """
        console_out = self.run_command('filtlong --target_bases 1 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_2'])
        self.assertTrue('target: 1 bp' in console_out)
        self.assertTrue('keeping 5,000 bp' in console_out)

    def test_sort_7(self):
        """
        Like test_sort_1, but with an assembly reference.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --target_bases 100000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 100,000 bp' in console_out)
        self.assertTrue('not enough reads to reach target' in console_out)

    def test_sort_8(self):
        """
        Like test_sort_2, but with an assembly reference.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --target_bases 10001 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 10,001 bp' in console_out)
        self.assertTrue('keeping 15,000 bp' in console_out)

    def test_sort_9(self):
        """
        Like test_sort_3, but with an assembly reference. With a reference, reads 1 and 3 are the best two.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --target_bases 10000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_3'])
        self.assertTrue('target: 10,000 bp' in console_out)
        self.assertTrue('keeping 10,000 bp' in console_out)

    def test_sort_10(self):
        """
        Like test_sort_3, but with an assembly reference. With a reference, reads 1 and 3 are the best two.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --target_bases 5001 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_3'])
        self.assertTrue('target: 5,001 bp' in console_out)
        self.assertTrue('keeping 10,000 bp' in console_out)

    def test_sort_11(self):
        """
        Like test_sort_5, but with an assembly reference. With a reference (instead of Phred), read 1 is the best.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --target_bases 5000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1'])
        self.assertTrue('target: 5,000 bp' in console_out)
        self.assertTrue('keeping 5,000 bp' in console_out)

    def test_sort_12(self):
        """
        Like test_sort_5, but with an assembly reference. With a reference (instead of Phred), read 1 is the best.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --target_bases 1 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1'])
        self.assertTrue('target: 1 bp' in console_out)
        self.assertTrue('keeping 5,000 bp' in console_out)

    def test_sort_13(self):
        """
        Like test_sort_7, but with a read reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --target_bases 100000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 100,000 bp' in console_out)
        self.assertTrue('not enough reads to reach target' in console_out)

    def test_sort_14(self):
        """
        Like test_sort_8, but with a read reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --target_bases 10001 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_2', 'test_sort_3'])
        self.assertTrue('target: 10,001 bp' in console_out)
        self.assertTrue('keeping 15,000 bp' in console_out)

    def test_sort_15(self):
        """
        Like test_sort_9, but with a read reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --target_bases 10000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_3'])
        self.assertTrue('target: 10,000 bp' in console_out)
        self.assertTrue('keeping 10,000 bp' in console_out)

    def test_sort_16(self):
        """
        Like test_sort_10, but with a read reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --target_bases 5001 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1', 'test_sort_3'])
        self.assertTrue('target: 5,001 bp' in console_out)
        self.assertTrue('keeping 10,000 bp' in console_out)

    def test_sort_17(self):
        """
        Like test_sort_11, but with a read reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --target_bases 5000 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1'])
        self.assertTrue('target: 5,000 bp' in console_out)
        self.assertTrue('keeping 5,000 bp' in console_out)

    def test_sort_18(self):
        """
        Like test_sort_12, but with a read reference.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --target_bases 1 INPUT > OUTPUT.fastq')
        output_reads = load_fastq(self.output_file)
        read_names = [x[0].decode() for x in output_reads]
        self.assertEqual(read_names, ['test_sort_1'])
        self.assertTrue('target: 1 bp' in console_out)
        self.assertTrue('keeping 5,000 bp' in console_out)
