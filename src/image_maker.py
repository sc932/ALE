#!/usr/bin/python

# (C) 2011 Scott Clark

import sys
import os
import subprocess
import commands
import time

error_script = './artificial_errors.py'
plotter = './plotter3.py'
original_file = 'Ecoli_first350k'

def run_it_through(file_name, error_opts):
    """Perform the following (if needed):

            1. make a fasta file using artificial_errors
            2. make a bowtie db
            3. run bowtie on the db and fasta file
            4. run ALE on the map and fasta
            5. run plotter3 on the ALE file
       """
    # make fasta file
    if not os.path.exists("%s.fna" % file_name):
        fasta_script = "%s %s %s.fna" % (error_script, error_opts, original_file)
        print fasta_script
        print commands.getoutput(fasta_script)
    # run bowtie on it
    if not os.path.exists("%s.1.ebwt" % file_name):
        bowtie_build = "bowtie-build %s.fna %s" % (file_name, file_name)
        print bowtie_build
        print commands.getoutput(bowtie_build)
    if not os.path.exists("%s.map.sam" % file_name):
        bowtie_script = "bowtie -t -I 0 -X 300 --fr -a -l 10 -v 1 -e 300 -S --threads 2 %s -1 part1_Ecoli_first350k.fastq  -2 part2_Ecoli_first350k.fastq %s.map.sam" % (file_name, file_name)
        print bowtie_script
        print commands.getoutput(bowtie_script)
    # run ALE on it
    if not os.path.exists("%s.ale" % file_name):
        ALE_script = "ALE %s.map.sam %s.fna %s.ale" % (file_name, file_name, file_name)
        print ALE_script
        print commands.getoutput(ALE_script)
    # run the plotter
    if not os.path.exists("%s.ale.pdf" % file_name):
        plot_script = "%s %s.ale" % (plotter, file_name)
        print plot_script
        print commands.getoutput(plot_script)

def main():
    """Perform the following (if needed):
            
            1. Download an E.Coli genome
            2. Truncate it to 350k bases
            3. Synthesize 2M reads
            4. Run ALE on some basic transformations
            5. Generate output plots
    """

    if not os.path.exists("CP000948.fna"):
        # wget from ncbi
        wget_command = "wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid20079/CP000948.fna"
        print wget_command
        print commands.getoutput(wget_command)
        
    if not os.path.exists("%s.fna" % original_file):
        # wget from ncbi and cat first 350k reads
        head_command = "head -50001 CP000948.fna > %s.fna" % original_file
        print head_command
        print commands.getoutput(head_command)

    if not os.path.exists("part1_%s.fastq" % original_file) or not os.path.exists("part2_%s.fastq" % original_file):
        # make reads
        synth_command = "./synthReadGen -ip 1.0 -nr 2000000 -ps 10 -b %s.fna %s.fastq" % (original_file, original_file)
        print synth_command
        print commands.getoutput(synth_command)

    ### sub/indel errors
    file_name = '%s_sub_errors' % original_file
    error_opts = "-ase 75000 2 -ade 175000 2 -aie 275000 2 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)
    
    ### transpose errors
    file_name = '%s_transpose_error' % original_file
    error_opts = "-trp 175000 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### inversion errors
    file_name = '%s_inversion_error' % original_file
    error_opts = "-inv 175000 200 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### sub/indel errors
    #file_name = '%s_broken_apart' % original_file
    #error_opts = "-ab 50000 -ab 150000 -o %s.fna" % file_name
    #run_it_through(file_name, error_opts)

    ### copy error
    file_name = '%s_copy_error' % original_file
    error_opts = "-cip 150000 50000 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

if __name__ == '__main__':
    main()
