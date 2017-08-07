'''
This script was a quick hack to make simulated Illumina and PacBio reads for circular sequences.
PBSIM and Wgsim do not seem to have a 'circular reference' option, nor did it have an option for
different depths for different references.  This script simply runs PBSIM over and over, each time
randomly rotating the reference sequences, adding all results to one file.

It's pretty crude - doesn't take command line arguments but instead you need to set constants at
the top of the script.
'''
from __future__ import print_function
from __future__ import division
import random
import os
import subprocess

# INPUT AND PARAMETERS
source_fasta = 'test_reference.fasta'
base_depths = [1.0] # One for each sequence in the source
pacbio_depth_factor = 30.0
illumina_depth_factor = 40.0
start_position_count = 10 # The number of times to run PBSIM and Wgsim with random start positions
illumina_read_length = 100
illumina_fragment_size = 350
illumina_fragment_std_dev = 50
pbsim_command = '/Users/Ryan/Applications/pbsim/pbsim' 
wgsim_command = 'wgsim'

# OUTPUT
pacbio_reads_filename = 'synthetic_pacbio_reads.fastq'
illumina_reads_filename_1 = 'synthetic_illumina_reads_1.fastq'
illumina_reads_filename_2 = 'synthetic_illumina_reads_2.fastq'

# Global var:
pacbio_read_counter = 1

def main():
    names, ref_seqs = load_fasta(source_fasta)
    pacbio_reads = open(pacbio_reads_filename, 'w')
    illumina_reads_1 = open(illumina_reads_filename_1, 'w')
    illumina_reads_2 = open(illumina_reads_filename_2, 'w')
    for i in range(start_position_count):
        for i, _ in enumerate(ref_seqs):
            generate_reads(ref_seqs[i], names[i], base_depths[i] / start_position_count,
                           pacbio_reads, illumina_reads_1, illumina_reads_2)
    pacbio_reads.close()
    illumina_reads_1.close()
    illumina_reads_2.close()
    zip_file(illumina_reads_filename_1)
    zip_file(illumina_reads_filename_2)
    zip_file(pacbio_reads_filename)

def generate_reads(sequence, name, base_depth, pacbio_reads, illumina_reads_1, illumina_reads_2):

    global pacbio_read_counter

    # Randomly rotate the sequence and save it to file.
    random_start = random.randint(0, len(sequence) - 1)
    rotated = sequence[random_start:] + sequence[:random_start]
    rotated_file = open('rotated.fasta', 'w')
    rotated_file.write('>' + name + '\n')
    rotated_file.write(rotated)
    rotated_file.write('\n')

    # Generate PacBio reads with PBSIM
    pacbio_depth = base_depth * pacbio_depth_factor
    command = [pbsim_command,
               '--data-typ', 'CLR',
               '--depth', str(pacbio_depth),
               '--model_qc', '/Users/Ryan/Applications/pbsim/model_qc_clr',
               'rotated.fasta']
    print(' '.join(command))
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    print(out)
    print(err)
    pbsim_file = open('sd_0001.fastq', 'r')

    # Copy reads to main file, but add the cycle number to the read ID so we don't end up with
    # multiple reads with the same ID
    line_num = 0
    for line in pbsim_file:
        if line_num % 2 == 0:
            line = line.strip() + '_' + str(pacbio_read_counter) + '\n'
        pacbio_reads.write(line)
        line_num += 1
        if line_num % 4 == 0:
            pacbio_read_counter += 1

    # Generate Illumina reads with Wgsim
    illumina_depth = base_depth * illumina_depth_factor
    illumina_read_pair_count = int(illumina_depth * len(sequence) / (2 * illumina_read_length))
    command = [wgsim_command,
               '-1', str(illumina_read_length),
               '-2', str(illumina_read_length),
               '-d', str(illumina_fragment_size),
               '-s', str(illumina_fragment_std_dev),
               '-N', str(illumina_read_pair_count),
               'rotated.fasta',
               '1.fastq',
               '2.fastq']
    print(' '.join(command))
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    print(out)
    print(err)
    wgsim_file_1 = open('1.fastq', 'r')
    for line in wgsim_file_1:
        illumina_reads_1.write(line)
    wgsim_file_2 = open('2.fastq', 'r')
    for line in wgsim_file_2:
        illumina_reads_2.write(line)

    # Clean up
    os.remove('sd_0001.fastq')
    os.remove('sd_0001.maf')
    os.remove('sd_0001.ref')
    os.remove('1.fastq')
    os.remove('2.fastq')
    os.remove('rotated.fasta')

def zip_file(filename):
    command = ['gzip', filename]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = process.communicate()
    
def load_fasta(filename): # type: (str) -> list[tuple[str, str]]
    '''
    Returns the names and sequences for the given fasta file.
    '''
    seq_names = []
    seqs = []
    fasta_file = open(filename, 'r')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>': # Header line = start of new contig
            if name:
                seq_names.append(name.split()[0])
                seqs.append(sequence)
                name = ''
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        seq_names.append(name.split()[0])
        seqs.append(sequence)
    return seq_names, seqs

if __name__ == '__main__':
    main()
