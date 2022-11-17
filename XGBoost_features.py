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
from operator import itemgetter
import numpy as np
from overlaps import parse_paf, filter_primary_overlaps, remove_overlap, Overlap, extend_overlaps
from sequences import *
#from data_parsers import *

from typing import *

#COV = int(1.2 * 35)
COV = 30
get_bases_in_order = itemgetter('A', 'C', 'G', 'T', 'D')

encode = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3
}

GAP_CODE = 4
def debug_overlaps(paf_path:str, reads_path:str, truth_path:str):
    print('Parsing inputs...')
    time1 = time.time()
    overlaps = parse_paf(paf_path)
    overlaps = filter_primary_overlaps(overlaps)
    extend_overlaps(overlaps)
    # overlaps = take_longest(overlaps)
    time2 = time.time()
    print(f'Time taken for parsing paf: {time2-time1}')
    covs = [len(olps) for olps in overlaps.values()]
    covs.sort()
    print('Median number of overlaps:', covs[len(covs) // 2])
    
    
    
    global reads, truth_reads
    time3 = time.time()
    reads = get_reads(reads_path)
    truth_reads = get_reads(truth_path)
    time4 = time.time()
    print(f"Time taken for parsing reads: {time4-time3}")
    return overlaps
def calculate_iden(cigar):
    matches, mis, ins, dels = 0, 0, 0, 0
    for op, length in cigar:
        if op == '=':
            matches += length
        elif op == 'X':
            mis += length
        elif op == 'I':
            ins += length
        elif op == 'D':
            dels += length
        else:
            ValueError('Invalid CIGAR')

    return matches / (matches + mis)  # + ins + dels)


class aux_data:
    def __init__(self, tname:str, ins_counts:List[int], num_pos:int):
        self.tname = tname
        
        # number of insertion slots at each position
        # of this target read
        self.ins_counts = ins_counts
        
        # total number of pos, including ins positions, 
        # of the alignment to this target read
        self.num_pos = num_pos
        
        
    def __repr__(self):
        return f'{self.tname} {self.num_pos}'


def one_hot(i):
    o = [0] * 5
    o[i] = 1
    return o

# returns a list with the original bases in each position (including ins pos) encoded 
def get_original_bases(original:str, ins_counts:List[int])->List[int]:
    list_of_lists=  [[encode[original[i]]] + count * [GAP_CODE] for i, count in enumerate(ins_counts)]
    return [item for sublist in list_of_lists for item in sublist]
'''
Returns 
    1. a list of list of dictionaries, containing the counts of bases in each positions 
    2. aux_data for this target read
'''
def get_bases_freq(reads: Dict[str, HAECSeqRecord], tname: str,
                  tovlps: List[Overlap]) -> Tuple[List[List[Dict[str, int]]], aux_data, List[int]]:
    # uncorrected
    uncorrected = reads[tname].seq
    # reads before error correction
    # seq_lst.append(reads[target_read])
    # print(len(uncorrected))
    # Add original base into frequency
    freq = []
    for i, base in enumerate(uncorrected):
        freq.append([])
        freq[i].append(defaultdict(int))
        freq[i][0][base] += 1
        #freq[i][0]['O'] = encode[base]

    # keep track of number of ins at each original position in uncorrected
    ins_counts = [0] * len(freq) 
    
    # keep track of total number of positions including ins (number of dict)
    total_num_pos = len(freq)
    # TODO check usefulness
    # temp_lst = []
    if len(tovlps) >= COV:
        tovlps.sort(key=lambda o: calculate_iden(o.cigar), reverse=True)
        #tovlps.sort(key=lambda o: o.iden, reverse=True)
        temp_lst = tovlps[:COV]
    else:
        temp_lst = tovlps
    #if tname == 'e012f204-6a49-4e82-884e-8138929a86c9_1':
    #    print('Lengths:', len(tovlps), len(temp_lst))
    tovlps = temp_lst

    # freq of A C G T and D at each base
    for overlap in tovlps:
        # sr = SeqRecord(reads[query_name].seq)
        # sr.id = reads[query_name].id
        # sr.name = reads[query_name].name
        # sr.description = reads[query_name].description
        # seq_lst.append(sr)
        #cnt = 0
        # reverse complement
        # recalculate the query start and end in case of RC
        # rc_start = len - end - 1
        # rc_end = len - start
        target_pos = overlap.tstart
        query_seq_record = reads[overlap.qname]
        if overlap.strand == '+':
            query_read = query_seq_record.seq
            query_pos = overlap.qstart
        if overlap.strand == '-':
            query_read = query_seq_record.reverse_complement()
            query_pos = len(query_read) - overlap.qend
            # query_end = len(query_read) - overlap.qstart
        # print(f'target_start: {target_start}, target_end:{target_end}')
        # print(uncorrected[target_start:target_end])
        # print(f'query_start{query_start}, query_end{query_end}')
        # print(query_read[query_start:query_end])

        for operation, length in overlap.cigar:
            if operation == '=' or operation == 'X':
                for i, b in enumerate(query_read[query_pos:query_pos + length]):
                    freq[target_pos + i][0][b] += 1

                target_pos += length
                query_pos += length
            elif operation == 'D':
                for i in range(target_pos, target_pos + length):
                    freq[i][0]['D'] += 1

                target_pos += length
            elif operation == 'I':
                target_tmp_pos = target_pos - 1
                for i, b in enumerate(query_read[query_pos:query_pos + length],
                                      start=1):
                    if i == len(freq[target_tmp_pos]):
                        freq[target_tmp_pos].append(defaultdict(int))
                        #freq[target_tmp_pos][i]['O'] = GAP_CODE
                        ins_counts[target_tmp_pos] += 1
                        total_num_pos += 1
                    freq[target_tmp_pos][i][b] += 1

                query_pos += length
            else:
                raise ValueError(f'{operation} - Invalid CIGAR operation.')
        # assert pointers are at the end
        # assert target_pos == target_end and query_pos == query_end, f'{target_start}, {target_pos}, {target_end}, {query_start}, {query_pos}, {query_end}, {strand}'
        # assert target_pos == target_end and query_pos == query_end, f'{target_read}, {query_name}'
        # assert target_pos == target_end and query_pos == query_end, f'{path}'
    # generate consensus
    
    return freq, aux_data(tname, ins_counts, total_num_pos), get_original_bases(uncorrected, ins_counts)
    
    

def generate_cigar(overlaps: Dict[str, List[Overlap]],
                   reads: Dict[str, HAECSeqRecord]) -> None:
    cnt = 0
    for tname, tovlps in overlaps.items():
        # if cnt == 1: break
        # print(cnt, "#######################")
        cnt += 1
        # print(len(overlap_map[query_seq]))
        count = 0
        for overlap in tovlps:
            if (overlap.tend - overlap.tstart) / len(reads[tname].seq) < 0.05:
                continue

            cigar = calculate_path(overlap, reads[tname], reads[overlap.qname])
            #overlap.cigar = calculate_iden(cigar)
            if cigar is not None:
                overlap.cigar = cigar
                count += 1
            else:
                pass
                #print(target_seq, overlap)
            # reverse_cigar = (query_start, query_end, target_seq, target_start, target_end, reverse(path), strand)
            # cigar_list[query_name].append(reverse_cigar)
        # if tname == 'e012f204-6a49-4e82-884e-8138929a86c9_1':
        #    print('Aln counts:', len(tovlps), count)


def calculate_path(overlap: Overlap, trecord: HAECSeqRecord,
                   qrecord: HAECSeqRecord):
    target_read = trecord.seq[overlap.tstart:overlap.tend]
    if overlap.strand == '-':
        query_read = qrecord.reverse_complement(overlap.qstart, overlap.qend)
    else:
        query_read = qrecord.seq[overlap.qstart:overlap.qend]

    align = edlib.align(query_read, target_read, task='path')
    path = align['cigar']
    # distance = align['editDistance']
    # print(distance/(target_end - target_start) * 100, '%')
    if path is None:
        # print(target_name, query_name, query_read, target_read.lower())
        return None

    generator = gen(path)
    cigar = list(generator)
    # dels = sum(i for _, i in path if _ == 'D')
    # inserts = sum(i for _, i in path if _ == 'I')
    # total_dels += dels
    # print(f'deletions: {dels}')
    # print(f'insertions: {inserts}')
    return cigar


#REVERSED_OP = {'D': 'I', 'I': 'D', '=': '=', 'X': 'X'}

#def reverse(path):
#    return [(REVERSED_OP[op], l) for (op, l) in reversed(path)]


def fix_gap_count(dicts_for_pos:List[Dict[str, int]])->None:
    n_support = sum(dicts_for_pos[0].values())
    for i in range(1, len(dicts_for_pos)):
        dicts_for_pos[i]['D'] = n_support - sum(dicts_for_pos[i].values()) #+ dicts_for_pos[i]['O']

def flatten_freq(freq:List[List[Dict[str, int]]])->List[Dict[str, int]]:
    for dicts_for_pos in freq:
        fix_gap_count(dicts_for_pos)
    return [item for sublist in freq for item in sublist]
    #flattened_freq = [item for sublist in freq for item in sublist]
    #return np.asarray([get_bases_in_order(counts_dict) for counts_dict in flattened_freq])
    

'''
This should return  
1. shape (n, 5) numpy arrays containing all
the positions of all the inputs sequences concatenated.
2. a list of aux_data, containing info on number of ins in each position of each target read
in the batch (for generating ground-truth) and total number of positions(dictionaries) for 
each read (for constructing corrected read from model output)
'''
def generate_input_features(reads: Dict[str, HAECSeqRecord], overlaps: Dict[str, List[Overlap]]):

    generate_cigar(overlaps, reads)

    freq_list = []
    aux_data_list = []
    all_original_bases = []
    for tname, tovlps in overlaps.items():
        tovlps = [o for o in tovlps if o.cigar is not None]
        freq, aux_data, original_bases = get_bases_freq(reads, tname, tovlps)
        aux_data_list.append(aux_data)
        freq_list += freq
        all_original_bases += original_bases
    
    # from list of list of dict to list of dict
    # also fixes gap count for ins positions
    flattened_freq = flatten_freq(freq_list)
    assert(len(all_original_bases) == len(flattened_freq))
    return np.asarray([list(get_bases_in_order(counts_dict)) + one_hot(original_base) for counts_dict, original_base in zip(flattened_freq, all_original_bases)]), aux_data_list

PATTERN = re.compile('(\d+)([=XID])')


def gen(string):
    for match in PATTERN.finditer(string):
        yield match.group(2), int(match.group(1))

 
'''
assumes the truth read is the same relative strand and has same name as uncorrected
'''
def ground_truth_for_read(aux : aux_data):
    global reads, truth_reads
    target_read = reads[aux.tname].seq
    truth_read = truth_reads[aux.tname].seq
    
    align = edlib.align(truth_read, target_read, task='path')
    path = align['cigar']
    
    ###
    '''
    nice = edlib.getNiceAlignment(align, truth_read, target_read)
    tr_align = nice['query_aligned']
    tg_align = nice['target_aligned']
    m_align = nice['matched_aligned']
    
    for i in range(0, len(tr_align), 80):
        print(tg_align[i:i+80])
        print(m_align[i:i+80])
        print(tr_align[i:i+80])
    '''
    
    ###
    
    generator = gen(path)
    cigar = list(generator)
    
    # list of label bases at each pos, encoded as 0,1,2,3,4
    truth_list = []
    truth_idx = 0
    uncorrected_idx = 0

    assert(cigar[0][0] != 'I')
    for op, length in cigar:
        if op == '=' or op == 'X':
            list_of_lists = [[encode[base]] + [GAP_CODE] * aux.ins_counts[uncorrected_idx+i] for i, base in enumerate(truth_read[truth_idx:truth_idx+length])]
            flattened_truth_values = [value for sublist in list_of_lists for value in sublist]
            truth_list += flattened_truth_values
            #truth_list += [encode[base] for base in truth_read[truth_idx:truth_idx+length]]
            truth_idx += length
            uncorrected_idx += length
            
        elif op == 'D':
            temp_list = [GAP_CODE] * (length + sum(aux.ins_counts[uncorrected_idx:uncorrected_idx+length]))
            truth_list += temp_list
            uncorrected_idx += length
            
            
        elif op == 'I':
            num_ins_pos_with_support = aux.ins_counts[uncorrected_idx-1]
            #print(f'unc {uncorrected_idx}, num_ins_pos_with: {num_ins_pos_with_support}')
            if num_ins_pos_with_support == 0:
                truth_idx += length
                continue
            temp_list = [encode[base] for base in truth_read[truth_idx:truth_idx+min(length, num_ins_pos_with_support)]]
            
            '''if num_ins_pos_with_support > length:
                temp_list += [GAP_CODE] * (num_ins_pos_with_support - length)
            truth_list[-num_ins_pos_with_support:] = temp_list[:]
            '''
            if num_ins_pos_with_support == len(temp_list):
                truth_list[-num_ins_pos_with_support:] = temp_list[:]
            else:
                truth_list[-num_ins_pos_with_support:-num_ins_pos_with_support+len(temp_list)] = temp_list[:]             
            
            truth_idx += length
            
            
        else:
            raise ValueError(f'{op} - Invalid CIGAR operation.')
    return truth_list       
'''
batch_input_features: (n, 5) np array
truth: (n, ) np array
'''
def generate_training_features(overlaps: Dict[str, List[Overlap]]):
    global reads
    batch_input_features, aux_data_list = generate_input_features(reads, overlaps)
    list_of_lists = [ground_truth_for_read(aux) for aux in aux_data_list]
    batch_ground_truth = [label for sub_list in list_of_lists for label in sub_list]
    return batch_input_features, np.asarray(batch_ground_truth)

'''
batch_input_features: (n, 5) np array
truth: (n, ) np array
'''
def debug_generate_training_features(overlaps: Dict[str, List[Overlap]]):
    global reads
    batch_input_features, aux_data_list = generate_input_features(reads, overlaps)
    list_of_lists = [ground_truth_for_read(aux) for aux in aux_data_list]
    batch_ground_truth = [label for sub_list in list_of_lists for label in sub_list]
    #read = reads['60b405f2-c66f-4cc4-8524-62e19b211c7f_1']
    return batch_input_features, np.asarray(batch_ground_truth), aux_data_list#, read

def take_longest(
        overlaps: Dict[str, List[Overlap]]) -> Dict[str, List[Overlap]]:
    longest = {}
    for tname, tovlps in overlaps.items():
        tovlps_to_qname = defaultdict(list)
        for o in tovlps:
            tovlps_to_qname[o.qname].append(o)

        for ovlps in tovlps_to_qname.values():
            ovlps.sort(key=lambda o: (o.tend - o.tstart), reverse=True)
        longest[tname] = [ovlps[0] for ovlps in tovlps_to_qname.values()]

    return longest



# Pool all the inputs features and ground-truth labels 
# for training
def generate_training_data(paf_path:str, reads_path:str, truth_path:str, num_proc:int, limit:int):
    print('Parsing inputs...')
    time1 = time.time()
    overlaps = parse_paf(paf_path)
    overlaps = filter_primary_overlaps(overlaps)
    extend_overlaps(overlaps)
    # overlaps = take_longest(overlaps)
    time2 = time.time()
    print(f'Time taken for parsing paf: {time2-time1}')
    covs = [len(olps) for olps in overlaps.values()]
    covs.sort()
    print('Median number of overlaps:', covs[len(covs) // 2])
    
    
    
    global reads, truth_reads
    time3 = time.time()
    reads = get_reads(reads_path)
    truth_reads = get_reads(truth_path)
    time4 = time.time()
    print(f"Time taken for parsing reads: {time4-time3}")
    

    futures_input_features = []
    overlap_keys = list(overlaps)
    overlap_list = overlaps.items()
    
    batch_input_features_list, batch_ground_truth_list = [], []
    time5 = time.time()
    if limit == -1:
        limit = len(overlap_list)
    with ProcessPoolExecutor(max_workers=num_proc) as executor:
        # step = int(len(overlap_map) / workers)
        step = 10
        # managed_reads = manager.dict(reads)
        # for i in range(0, 1000, step):
        for i in range(0, limit, step):
            end = min(i + step, limit)
            curr_dict = {k: overlaps[k] for k in overlap_keys[i:end]}

            f = executor.submit(generate_training_features, curr_dict)
            futures_input_features.append(f)

        for result in tqdm(as_completed(futures_input_features)):
            # seq_lst.extend(result.result())
            batch_input_features, batch_ground_truth = result.result()
            batch_input_features_list.append(batch_input_features)
            batch_ground_truth_list.append(batch_ground_truth)
            
    all_input_features = np.concatenate(batch_input_features_list)
    all_ground_truth = np.concatenate(batch_ground_truth_list)
    time6 = time.time()
    print(f'Time taken to generate training features for {limit} targets: {time6-time5}')
    return all_input_features, all_ground_truth
        
