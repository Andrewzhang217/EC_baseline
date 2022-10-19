#!/usr/bin/env python
from asyncio import as_completed
import multiprocessing
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from tqdm import tqdm
from pathlib import Path
from collections import defaultdict, Counter
import time
import argparse
from multiprocessing import set_start_method, Manager
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
import itertools
import re
import sys
from functools import wraps  
    

class CustomSeqRecord:
    RC_BASE = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    def __init__(self, name, id, description, seq):
        self.name = name
        self.id = id
        self.description = description
        self.seq = seq
        
    def reverse_complement(self):
        result = ''.join([(self.RC_BASE[base]) for base in reversed(self.seq)])
        return result


def correct_error(reads, target_read, lstDict):

    # uncorrected
    uncorrected = reads[target_read].seq
    # reads before error correction
    # seq_lst.append(reads[target_read])
    # print(len(uncorrected))

    # Add original base into frequency
    freq = []
    for i, base in enumerate(uncorrected):
        freq.append([])
        freq[i].append(defaultdict(int))
        freq[i][0][base] += 1


    # freq of A C G T and D at each base
    for target_start, target_end, query_name, query_start, query_end, path, strand in lstDict:
        # sr = SeqRecord(reads[query_name].seq)
        # sr.id = reads[query_name].id
        # sr.name = reads[query_name].name
        # sr.description = reads[query_name].description
        # seq_lst.append(sr)
        cnt = 0
        # reverse complement
        # recalculate the query start and end in case of RC
        # rc_start = len - end - 1
        # rc_end = len - start
        target_pos = target_start
        query_pos = query_start
        query_seq_record = reads[query_name]
        query_read = query_seq_record.seq
        if strand == '-':
            query_read = query_seq_record.reverse_complement()
            query_pos = len(query_read) - query_end
            query_end = len(query_read) - query_start
        # print(f'target_start: {target_start}, target_end:{target_end}')
        # print(uncorrected[target_start:target_end])
        # print(f'query_start{query_start}, query_end{query_end}')
        # print(query_read[query_start:query_end])
        
        for operation, length in path:
            if operation == '=' or operation == 'X':
                for i, b in enumerate(query_read[query_pos:query_pos + length]):
                    freq[target_pos+i][0][b] += 1
                
                target_pos += length
                query_pos += length
            elif operation == 'D':
                for i in range(target_pos, target_pos + length):
                    freq[i][0]['D'] += 1

                target_pos += length    
            elif operation == 'I':
                target_tmp_pos = target_pos - 1
                for i, b in enumerate(query_read[query_pos:query_pos + length], start=1):
                    if i == len(freq[target_tmp_pos]):
                        freq[target_tmp_pos].append(defaultdict(int))
                    freq[target_tmp_pos][i][b] += 1

                query_pos += length
            else:
                raise ValueError(f'{operation} - Invalid CIGAR operation.')  
        # assert pointers are at the end
        assert target_pos == target_end and query_pos == query_end, f'{target_start}, {target_pos}, {target_end}, {query_start}, {query_pos}, {query_end}, {strand}'
        # assert target_pos == target_end and query_pos == query_end, f'{target_read}, {query_name}'
        # assert target_pos == target_end and query_pos == query_end, f'{path}'
    # generate consensus
    corrected = generate_consensus(uncorrected, freq, target_read)
    return corrected


def generate_consensus(uncorrected, freq, r_id):
    corrected = ''
    # counter = 0
    # dels = 0
    # haps = 0
    insertions, dels = 0, 0
    for pos in range(len(freq)):
        n_support = sum(freq[pos][0].values())

        # if r_id == '6c2efce6-9b99-4126-8444-0f32e8366b1d' and 955 <= len(corrected) < 966:
        #     print(freq[pos]) 
        for i in range(len(freq[pos])):
            if i > 0:
                freq[pos][i]['D'] = n_support - sum(freq[pos][i].values())

            # max_occ = max(freq[pos][i].values(), default=0)
            mc = (Counter(freq[pos][i]).most_common(2))
            if len(mc) == 1:
                b1, c1 = mc[0]
                c2 = 0
            else:
                (b1, c1), (_, c2) = mc[0], mc[1]
            # print(max_occ)
            if c1 == 0:
                break
            if c1 < 5:
                # #print(cnt)
                # counter += 1
                if i == 0:
                    corrected += uncorrected[pos]
            else :  # bimodal
                if c1 <= 1.5 * c2:
                    if i == 0:
                        corrected += uncorrected[pos]
                    else:
                        if b1 != 'D':
                            corrected += b1
                else:  # Unimodal
                    if b1 != 'D':
                        if i != 0:
                            insertions += 1
                        corrected += b1
                    else:
                        if i == 0:
                            dels += 1

    # print("insert" , insertions, dels)                    
    return corrected


def generate_cigar(overlap_map, reads):
    cigar_list = defaultdict(list)
    cnt = 0
    for target_name in overlap_map: 
        # if cnt == 1: break
        # print(cnt, "#######################")
        cnt += 1
        # print(len(overlap_map[query_seq]))
        for overlap in overlap_map[target_name]:
            aligned_overlap = calculate_path(overlap, reads, target_name)
            if aligned_overlap is not None:
                cigar_list[target_name].append(aligned_overlap)
            else:
                pass
                #print(target_seq, overlap)
            # reverse_cigar = (query_start, query_end, target_seq, target_start, target_end, reverse(path), strand)
            # cigar_list[query_name].append(reverse_cigar)

    return cigar_list


def calculate_path(overlap, reads, target_name):
    target_start = overlap[0]
    target_end = overlap[1]
    target_read = reads[target_name].seq[target_start:target_end]
    strand = overlap[2]
    query_name = overlap[3]
    query_start = overlap[4]
    query_end = overlap[5]
    if strand == '-':
        query_read = reads[query_name].reverse_complement()[query_start:query_end]
    else :
        query_read = reads[query_name].seq[query_start:query_end]
    align = edlib.align(query_read, target_read, task = 'path')
    path = align['cigar']
    # distance = align['editDistance']
    # print(distance/(target_end - target_start) * 100, '%')
    if path is None:
        print(target_name, query_name, query_read, target_read.lower())
        return None

    generator = gen(path)
    path = list(generator)
    # dels = sum(i for _, i in path if _ == 'D')
    # inserts = sum(i for _, i in path if _ == 'I')
    # total_dels += dels
    # print(f'deletions: {dels}')
    # print(f'insertions: {inserts}')
    return target_start, target_end, query_name, query_start, query_end, path, strand


REVERSED_OP = {'D': 'I', 'I': 'D', '=': '=', 'X': 'X'}
def reverse(path):
    return [(REVERSED_OP[op], l) for (op, l) in reversed(path)]

def gen(string):
    pattern = re.compile('(\d+)([=XID])')
    for match in pattern.finditer(string):
        yield match.group(2), int(match.group(1))

def get_reads(input_file):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'
    
    records_map = {}
    # records = list(SeqIO.parse(input_file, file_type))
    for record in SeqIO.parse(input_file, file_type):
        record.seq = str(record.seq).upper()
        records_map[record.id] = CustomSeqRecord(record.name, record.id, record.description, record.seq)

    print("size: ", len(records_map))
    return records_map 


def parse_paf(paf_file):
    overlap_map = defaultdict(list)
    with open(paf_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            
            query_name = line[0]
            query_start = int(line[2])
            query_end = int(line[3])
            strand = line[4] 
            target_name = line[5]
            target_start = int(line[7])
            target_end = int(line[8])
        
            #   parse to correct data type here
            overlap_map[query_name].append((query_start, query_end, strand, target_name, target_start, target_end))
    
    return overlap_map


def parse_args():
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-p', '--paf', type=str, help='path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    parser.add_argument('-t', '--thread', type=int, help='number of threads(processes actually')
    args = parser.parse_args()
    return args


def generate_cigar_correct_error(overlap_map):
    global reads

    cigar_list = generate_cigar(overlap_map, reads)
    seq_lst = []
    for target_read in cigar_list:
        corrected = correct_error(reads, target_read, cigar_list[target_read])
        corrected_seq = SeqRecord(Seq(corrected))
        corrected_seq.id = reads[target_read].id
        corrected_seq.name = reads[target_read].name
        corrected_seq.description = reads[target_read].description
        seq_lst.append(corrected_seq) 

    return seq_lst          


def main(args):
    # reads = get_reads(args.input)

    overlap_map = parse_paf(args.paf)
    print("finish parsing")
    workers = args.thread
    seq_lst = []
    futures_ec_reads = []
    overlap_keys = list(overlap_map)
    overlap_list = overlap_map.items()
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        # step = int(len(overlap_map) / workers)
        step = 10
        # managed_reads = manager.dict(reads)
        # for i in range(0, 1000, step):
        for i in range(0, len(overlap_list), step):
            end = min(i + step, len(overlap_list))
            curr_dict = {k : overlap_map[k] for k in overlap_keys[i : end]}
            f = executor.submit(generate_cigar_correct_error, curr_dict)
            futures_ec_reads.append(f)

        for result in tqdm(as_completed(futures_ec_reads)):
            seq_lst.extend(result.result())
    SeqIO.write(seq_lst, args.output, 'fasta')


def wrapper():
    set_start_method('forkserver')
   
    args = parse_args()
    
    t1 = time.time()
    main(args)
    t2 = time.time()
    print('Time elapsed: ', t2 - t1)


if __name__ == '__main__':
    args = parse_args()
    reads = get_reads(args.input)
    # set_start_method('forkserver')
    
    t1 = time.time()
    main(args)
    t2 = time.time()
    print('Time elapsed: ', t2 - t1)
