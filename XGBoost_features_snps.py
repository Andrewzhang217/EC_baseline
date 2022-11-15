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
from align_info import read_records

from data_parsers_snps import *

from typing import *

COV = int(1.2 * 35)
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
    def __init__(self, tname:str, ins_counts:list[int], num_pos:int):
        self.tname = tname
        
        # number of insertion slots at each position
        # of this target read
        self.ins_counts = ins_counts
        
        # total number of pos, including ins positions, 
        # of the alignment to this target read
        self.num_pos = num_pos

        
    def __repr__(self):
        return f'{self.tname} {self.num_pos}'

'''
Returns 
    1. a list of list of dictionaries, containing the counts of bases in each positions 
    2. aux_data for this target read
'''
def one_hot(i):
    o = [0] * 5
    o[i] = 1
    return o


def count_freq_low_fp(reads: Dict[str, HAECSeqRecord], tovlps: List[Overlap],
                      target_seq_record):
        # freq of A C G T and D at each base
    freq = []
    uncorrected = target_seq_record.seq
    for i, base in enumerate(uncorrected):
        freq.append([])
        freq[i].append(defaultdict(int))
        freq[i][0][base] += 1
    aux = aux_data(target_seq_record.id, [0] * len(freq), 0)
    for overlap in tovlps:
    
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
                        aux.ins_counts[target_tmp_pos] += 1
                        aux.num_pos += 1
                    freq[target_tmp_pos][i][b] += 1

                query_pos += length
            else:
                raise ValueError(f'{operation} - Invalid CIGAR operation.')
        return freq, aux, [[defaultdict(int)] * len(sublist) for sublist in freq]
def count_freq_high_fp(reads: Dict[str, HAECSeqRecord], tovlps: List[Overlap],
                      target_seq_record):
            # freq of A C G T and D at each base
    freq = []
    uncorrected = target_seq_record.seq
    for i, base in enumerate(uncorrected):
        freq.append([])
        freq[i].append(defaultdict(int))
        freq[i][0][base] += 1
    fp_freq = [[defaultdict(int)] for ls in freq]
    aux = aux_data(target_seq_record.id, [0] * len(freq), 0)
    for overlap in tovlps:
        
        target_pos = overlap.tstart
        query_seq_record = reads[overlap.qname]
        is_fp = target_seq_record.info.overlaps_with(query_seq_record.info)
        if overlap.strand == '+':
            query_read = query_seq_record.seq
            query_pos = overlap.qstart
        if overlap.strand == '-':
            query_read = query_seq_record.reverse_complement()
            query_pos = len(query_read) - overlap.qend

        for operation, length in overlap.cigar:
            if operation == '=' or operation == 'X':
                for i, b in enumerate(query_read[query_pos:query_pos + length]):
                    freq[target_pos + i][0][b] += 1
                    if is_fp:
                        fp_freq[target_pos + i][0][b] += 1

                target_pos += length
                query_pos += length
            elif operation == 'D':
                for i in range(target_pos, target_pos + length):
                    freq[i][0]['D'] += 1
                    if is_fp:
                        fp_freq[i][0]['D'] += 1
                target_pos += length
            elif operation == 'I':
                target_tmp_pos = target_pos - 1
                for i, b in enumerate(query_read[query_pos:query_pos + length],
                                      start=1):
                    if i == len(freq[target_tmp_pos]):
                        freq[target_tmp_pos].append(defaultdict(int))
                        #freq[target_tmp_pos][i]['O'] = GAP_CODE
                        if is_fp:
                            fp_freq[target_tmp_pos].append(defaultdict(int))
                        aux.ins_counts[target_tmp_pos] += 1
                        aux.num_pos += 1
                    freq[target_tmp_pos][i][b] += 1

                query_pos += length
            else:
                raise ValueError(f'{operation} - Invalid CIGAR operation.')
        return freq, aux, fp_freq
def count_freq(reads: Dict[str, HAECSeqRecord], tovlps: List[Overlap],
                      target_seq_record):
        # freq of A C G T and D at each base
    freq = []
    uncorrected = target_seq_record.seq
    for i, base in enumerate(uncorrected):
        freq.append([])
        freq[i].append(defaultdict(int))
        freq[i][0][base] += 1
    aux = aux_data(target_seq_record.id, [0] * len(freq), 0)
    for overlap in tovlps:
    
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
                        aux.ins_counts[target_tmp_pos] += 1
                        aux.num_pos += 1
                    freq[target_tmp_pos][i][b] += 1

                query_pos += length
            else:
                raise ValueError(f'{operation} - Invalid CIGAR operation.')
        return freq, aux, None
def num_fp_gt(tovlps: List[Overlap], num:int, target_seq_record:HAECSeqRecord):
    count = 0
    for overlap in tovlps:
        query_seq_record = reads[overlap.qname]
        if not target_seq_record.info.overlaps_with(query_seq_record.info):
            count+=1
    return count > num
def get_bases_freq_snps(reads: Dict[str, HAECSeqRecord], tname: str,
                  tovlps: List[Overlap]) -> list[list[dict[str, int]]]:
    # uncorrected
    target_seq_record = reads[tname]
    # TODO check usefulness
    # temp_lst = []
    if len(tovlps) >= COV:
        tovlps.sort(key=lambda o: calculate_iden(o.cigar), reverse=True)
        #tovlps.sort(key=lambda o: o.iden, reverse=True)
        temp_lst = tovlps[:COV]
    else:
        temp_lst = tovlps
    tovlps = temp_lst

    '''
    TODO
    keep a count of 
    num of fp overlaps. Go through overlaps once to count. 
    If number > some k, when doing the iteration below, keep for each position another freq table
    with the counts for the fp overlapping reads.
    
    Later will use this freq and the overall freq to determine which pos have snp.
    ''' 
    fp_freq = None
    if target_seq_record.info:
        fp_gt = num_fp_gt(tovlps, 5, target_seq_record)
        if fp_gt:
            freq, aux, fp_freq = count_freq_high_fp(reads, tovlps, target_seq_record)
        else:
            freq, aux, fp_freq = count_freq_low_fp(reads, tovlps, target_seq_record)
    else:
        freq, aux, fp_freq = count_freq(reads, tovlps, target_seq_record)       
        
        # assert pointers are at the end
        # assert target_pos == target_end and query_pos == query_end, f'{target_start}, {target_pos}, {target_end}, {query_start}, {query_pos}, {query_end}, {strand}'
        # assert target_pos == target_end and query_pos == query_end, f'{target_read}, {query_name}'
        # assert target_pos == target_end and query_pos == query_end, f'{path}'
    # generate consensus
    return freq, aux, fp_freq
    

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


def fix_gap_count(dicts_for_pos):
    n_support = sum(dicts_for_pos[0].values())
    for i in range(1, len(dicts_for_pos)):
        dicts_for_pos[i]['D'] = n_support - sum(dicts_for_pos[i].values()) #+ dicts_for_pos[i]['O']

def flatten_freq(freq):
    for dicts_for_pos in freq:
        fix_gap_count(dicts_for_pos)
    return [item for sublist in freq for item in sublist]
    #flattened_freq = [item for sublist in freq for item in sublist]
    #return np.asarray([get_bases_in_order(counts_dict) for counts_dict in flattened_freq])
    
def get_snp_GT_list(freq_number_lists, fp_freq_number_lists):
    return None

'''
This should return  
1. shape (n, 5) numpy arrays containing all
the positions of all the inputs sequences concatenated.
2. a list of aux_data, containing info on number of ins in each position of each target read
in the batch (for generating ground-truth) and total number of positions(dictionaries) for 
each read (for constructing corrected read from model output)
'''
def generate_input_features(reads, overlaps: Dict[str, List[Overlap]]):
    generate_cigar(overlaps, reads)

    freq_list = []
    fp_freq_list = []
    aux_data_list = []
    for tname, tovlps in overlaps.items():
        tovlps = [o for o in tovlps if o.cigar is not None]
        freq, aux_data, fp_freq = get_bases_freq_snps(reads, tname, tovlps)
        aux_data_list.append(aux_data)
        freq_list += freq
        fp_freq_list += fp_freq
    '''
    TODO 
    the fp_freqs may be None. in which case there's no truth 
    '''
    
    # from list of list of dict to list of dict
    # also fixes gap count for ins positions
    flattened_freq = flatten_freq(freq_list)
    flattened_fp_freq = flatten_freq(fp_freq_list)
    freq_number_lists = [get_bases_in_order(counts_dict) for counts_dict in flattened_freq]
    fp_freq_number_lists = [get_bases_in_order(counts_dict) for counts_dict in flattened_fp_freq]
    
    return np.asarray(freq_number_lists), aux_data_list

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
    read = reads['60b405f2-c66f-4cc4-8524-62e19b211c7f_1']
    return batch_input_features, np.asarray(batch_ground_truth), aux_data_list, read

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
def generate_training_data(paf_path:str, reads_path:str, align_info_path:str, num_proc:int, limit:int):
    print('Parsing inputs...')
    time1 = time.time()
    overlaps = parse_paf(paf_path)
    # overlaps = take_longest(overlaps)
    time2 = time.time()
    print(f'Time taken for parsing paf: {time2-time1}')
    covs = [len(olps) for olps in overlaps.values()]
    covs.sort()
    print('Median number of overlaps:', covs[len(covs) // 2])
    
    
    
    global reads
    time3 = time.time()
    align_info_records = read_records(align_info_path)
    reads = get_reads_snps(reads_path, align_info_records)
    del align_info_records
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
        

def parse_args():
    parser = argparse.ArgumentParser(
        description='Feature generation in baseline model style for XGBoost.')
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='path of the input FASTA/Q file.',
                        required=True)
    parser.add_argument('-p',
                        '--paf',
                        type=str,
                        help='path of the input PAF file.',
                        required=True)
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        help='number of threads(processes actually',
                        default=1)
    args = parser.parse_args()
    return args

'''
def wrapper():
    set_start_method('forkserver')

    args = parse_args()

    t1 = time.time()
    main(args)
    t2 = time.time()
    print('Time elapsed: ', t2 - t1)
'''

def test(paf_path, input_path, truth_path):
    global reads, truth_reads
    truth_reads = get_reads(truth_path)
    reads = get_reads(input_path)
    overlaps = parse_paf(paf_path)
    # overlaps = take_longest(overlaps)

    print("finish parsing")
    covs = [len(olps) for olps in overlaps.values()]
    covs.sort()
    print('Median number of overlaps:', covs[len(covs) // 2])
    return overlaps

    

def test2(overlaps, i, step):
    overlap_keys = list(overlaps)
    overlap_list = overlaps.items()
    end = min(i + step, len(overlap_list))
    curr_dict = {k: overlaps[k] for k in overlap_keys[i:end]}
    time1 = time.time()
    out =  generate_training_features(curr_dict)
    time2 = time.time()
    print(f'time taken for 10: {time2-time1}')
    return out
'''
if __name__ == '__main__':
    args = parse_args()
    reads = get_reads(args.input)
    # set_start_method('forkserver')

    t1 = time.time()
    main(args)
    t2 = time.time()
    print('Time elapsed: ', t2 - t1)

'''
