#!/usr/bin/env python


import argparse
import os
from pathlib import Path
import pysam
from run_programs import run
from sys import stderr
from typing import *
'''
Inputs:
1. Two sets of reads from different haplotypes
2. Two assemblies for the haplotypes

Map each set of reads to its respectively haplotype's assembly.
For each mapped read, output a record containing:
1. read id
2. mapped target
3. strand
4. query_start
5. query_end
6. reference_start
7. reference_end
'''

class align_info:
    def __init__(self, tname:str, is_reverse:bool, q_start:int, q_end:int, t_start:int, t_end:int):
        self.tname = tname 
        self.is_reverse = is_reverse
        self.q_start = q_start
        self.q_end = q_end
        self.t_start = t_start
        self.t_end = t_end
    
    def overlaps_with(self, other):
        if self.tname != other.tname:
            return False
        return self.t_start <= other.t_end and other.t_start <= self.t_end
    
    def __repr__(self):
        if self.is_reverse:
            return f'{self.tname}\t-\t{self.q_start}\t{self.q_end}\t{self.t_start}\t{self.t_end}'
        return f'{self.tname}\t+\t{self.q_start}\t{self.q_end}\t{self.t_start}\t{self.t_end}'


class align_record(align_info):
    def __init__(self, qname:str, tname:str, is_reverse:bool, q_start:int, q_end:int, t_start:int, t_end:int):
        self.qname = qname
        super.__init__(tname, is_reverse, q_start, q_end, t_start, t_end)
        '''self.tname = tname 
        self.is_reverse = is_reverse
        self.q_start = q_start
        self.q_end = q_end
        self.t_start = t_start
        self.t_end = t_end'''
    def __repr__(self):
        if self.is_reverse:
            return f'{self.qname}\t{self.tname}\t-\t{self.q_start}\t{self.q_end}\t{self.t_start}\t{self.t_end}'
        return f'{self.qname}\t{self.tname}\t+\t{self.q_start}\t{self.q_end}\t{self.t_start}\t{self.t_end}'



def align_reads(mm2_path:str, ref_path:str, reads_path:str, output_sam_path:str, num_threads:int, reuse:bool):
    """
    Align the reads with minimap2, output the sam into the output path.
    """
    if os.path.isfile(output_sam_path):
        print(f'{output_sam_path} already exists!', file=stderr)
        if not reuse:
            exit(0)
        return
    
    named_arguments = {
        '-t' : str(num_threads)
    }
    
    single_arguments = ['-a', ref_path, reads_path]
    with open(output_sam_path, 'w') as output_sam_handle:
        run(mm2_path, named_arguments, single_arguments, stdout=output_sam_handle)



def extract_alignments(input_sam_path:str):
    """
    Given a sam, 
    for each record, output a record
    containing: query name, target name, strand, q_start, q_end, t_start, t_end
    """
    with pysam.AlignmentFile(input_sam_path, "r") as samfile:
        for record in samfile.fetch():
            # without index will use until_eof=True, so will have unmapped
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            
            target_name = record.reference_name
            target_align_start = record.reference_start
            target_align_end = record.reference_end
            
            query_name = record.query_name
            query_align_start = record.query_alignment_start
            query_align_end = record.query_alignment_end
            
            if record.is_reverse:
                query_len = record.infer_read_length()
                yield align_record(query_name, target_name, True, query_len - query_align_end, 
                                   query_len - query_align_start, target_align_start, target_align_end)
            else:
                yield align_record(query_name, target_name, False, query_align_start, query_align_end,
                                   target_align_start, target_align_end)

def output_records(records, output_handle):
    for record in records:
        output_handle.write(record.__repr__() + '\n')


def prepare_dir(dir_path:str, reuse:bool):
    if os.path.isdir(dir_path):
        print(f'{dir_path} already exists!', file=stderr)
        if not reuse:
            exit(0)
    else:
        os.mkdir(dir_path)
        
def process_line(line:str):
    line = line.strip().split('\t')
    if line[2] == '+':
        is_reverse = False
    elif line[2] == '-':
        is_reverse = True
    else:
        raise ValueError(f'{line[2]} is not a valid relative strand!')
    return line[0], align_info(line[1], is_reverse, line[3], line[4], line[5], line[6])
# read records and return a dictionary with the query read ids as key        
def read_records(record_file_path:str) -> Dict[str, align_record]:
    with open(record_file_path, 'r') as f:
        return {qname : info for qname, info in [process_line(line) for line in f]}

def parse_args():
    parser = argparse.ArgumentParser(description='Extract the alignments of reads to a reference.')
    parser.add_argument('--i1', type=str, required=True, help='Path to reads of haplotype 1.')
    parser.add_argument('--i2', type=str, required=True, help='Path to reads of haplotype 2.')
    parser.add_argument('--a1', type=str, required=True, help='Path to hap1 assembly.')
    parser.add_argument('--a2', type=str, required=True, help='Path to hap2 assembly.')
    parser.add_argument('-m', '--mm2', type=str, default='minimap2', help='Path to minimap2 executable, will search in PATH if not specified.')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory path.')
    parser.add_argument('--reuse', action='store_true', help='Reuse existing files.')
    return parser.parse_args()

'''
TODO filter the reads to only those with an alignment on the refs. 
So that every read in the training set have an 'origin' on the genome defined.
Also, maybe clip the reads to only the aligned portions...
--> actually the tagged clipped reads already have them filtered, and the filtering 
should not differ if the parameters etc are the same.
'''

def main():
    args = parse_args()
    
    prepare_dir(args.output, args.reuse)
    
    align_reads_sam_path1 = os.path.join(args.output, Path(args.i1).stem + '_to_' + Path(args.a1).stem + '.sam')
    align_reads_sam_path2 = os.path.join(args.output, Path(args.i2).stem + '_to_' + Path(args.a2).stem + '.sam')
    
    align_reads(args.mm2, args.a1, args.i1, align_reads_sam_path1, args.threads, args.reuse)
    align_reads(args.mm2, args.a2, args.i2, align_reads_sam_path2, args.threads, args.reuse)
    
    records_output_path = os.path.join(args.output, 'records')
    with open(records_output_path, 'w') as records_output_handle: 
        gen1 = extract_alignments(align_reads_sam_path1)
        output_records(gen1, records_output_handle)
        gen2 = extract_alignments(align_reads_sam_path2)
        output_records(gen2, records_output_handle)

if __name__ == '__main__':
    main()