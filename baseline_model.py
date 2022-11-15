#!/usr/bin/env python
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
from overlaps import parse_paf, filter_primary_overlaps, remove_overlap, Overlap, extend_overlaps
from sequences import *
#from data_parsers import *

from typing import *

COV = int(1.2 * 35)


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


def correct_error(reads: Dict[str, HAECSeqRecord], tname: str,
                  tovlps: List[Overlap]) -> str:
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

    # TODO check usefulness
    # temp_lst = []
    if len(tovlps) >= COV:
        tovlps.sort(key=lambda o: calculate_iden(o.cigar), reverse=True)
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
        cnt = 0
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
                    freq[target_tmp_pos][i][b] += 1

                query_pos += length
            else:
                raise ValueError(f'{operation} - Invalid CIGAR operation.')
        # assert pointers are at the end
        # assert target_pos == target_end and query_pos == query_end, f'{target_start}, {target_pos}, {target_end}, {query_start}, {query_pos}, {query_end}, {strand}'
        # assert target_pos == target_end and query_pos == query_end, f'{target_read}, {query_name}'
        # assert target_pos == target_end and query_pos == query_end, f'{path}'
    # generate consensus
    corrected = generate_consensus(uncorrected, freq)
    return corrected


def generate_consensus(uncorrected, freq):
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
            else:  # bimodal
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

PATTERN = re.compile('(\d+)([=XID])')


def gen(string):
    for match in PATTERN.finditer(string):
        yield match.group(2), int(match.group(1))


def generate_cigar_correct_error(overlaps: Dict[str, List[Overlap]]):
    global reads

    generate_cigar(overlaps, reads)

    seq_lst = []
    for tname, tovlps in overlaps.items():
        tovlps = [o for o in tovlps if o.cigar is not None]
        corrected = correct_error(reads, tname, tovlps)

        # TODO to HAECRecord
        corrected_seq = SeqRecord(Seq(corrected))
        corrected_seq.id = reads[tname].id
        corrected_seq.name = reads[tname].name
        corrected_seq.description = reads[tname].description
        seq_lst.append(corrected_seq)

    return seq_lst


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


def main(args):
    overlaps = parse_paf(args.paf)
    overlaps = filter_primary_overlaps(overlaps)
    extend_overlaps(overlaps)
    # overlaps = take_longest(overlaps)

    print("finish parsing")
    covs = [len(olps) for olps in overlaps.values()]
    covs.sort()
    print('Median number of overlaps:', covs[len(covs) // 2])

    futures_ec_reads = []
    overlap_keys = list(overlaps)
    overlap_list = overlaps.items()

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        # step = int(len(overlap_map) / workers)
        step = 10
        # managed_reads = manager.dict(reads)
        # for i in range(0, 1000, step):
        for i in range(0, len(overlap_list), step):
            end = min(i + step, len(overlap_list))
            curr_dict = {k: overlaps[k] for k in overlap_keys[i:end]}

            f = executor.submit(generate_cigar_correct_error, curr_dict)
            futures_ec_reads.append(f)

        with open(args.output, 'w') as f_out:
            for result in tqdm(as_completed(futures_ec_reads)):
                # seq_lst.extend(result.result())
                SeqIO.write(result.result(), f_out, 'fasta')

            # Write uncorrected reads
            for r_id, record in reads.items():
                if r_id not in overlaps:
                    f_out.write(f'>{r_id} {record.description}\n')
                    f_out.write(f'{record.seq}\n')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Baseline model for error correction')
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
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='path of error corrected reads file',
                        required=True)
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        help='number of threads(processes actually',
                        default=1)
    args = parser.parse_args()
    return args


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
