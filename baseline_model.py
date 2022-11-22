from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from tqdm import tqdm
import numpy as np

from collections import defaultdict, Counter
import time
import argparse
from multiprocessing import set_start_method
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
import sys

from sequences import *
from overlaps import parse_paf, filter_primary_overlaps, remove_overlap, Overlap, extend_overlaps, remove_overlap_ratio

from typing import *

# COV = int(1.2 * 35)
MIN_BASE = 5


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

    return matches / (matches + mis + ins + dels)


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
    # tovlps = [o for o in tovlps if calculate_iden(o.cigar) >= 0.9]

    # k = int(30 * len(reads[tname].seq) / median_ovlp_length)
    '''k = 30
    if len(tovlps) >= k:
        tovlps.sort(key=lambda o:
                    (o.tend - o.tstart) / o.tlen * calculate_iden(o.cigar)**2,
                    reverse=True)
        temp_lst = tovlps[:k]
    else:
        temp_lst = tovlps
    tovlps = temp_lst'''

    # freq of A C G T and D at each base
    ovlp_pos = set()
    for overlap in tovlps:
        target_pos = overlap.tstart
        query_seq_record = reads[overlap.qname]
        if overlap.strand == '+':
            query_read = query_seq_record.seq
            query_pos = overlap.qstart
            #query_end = overlap.qend
        if overlap.strand == '-':
            query_read = query_seq_record.reverse_complement()
            query_pos = len(query_read) - overlap.qend
            #query_end = len(
            #    query_read) - overlap.qstart  # TODO Comment after asserting

        for cig_pos, (operation, length) in enumerate(overlap.cigar):
            if operation == '=' or operation == 'X':
                for i, b in enumerate(query_read[query_pos:query_pos + length]):
                    freq[target_pos + i][0][b] += 1
                    if target_pos + i == 177637:
                        ovlp_pos.add((overlap.qname, query_pos + i))

                target_pos += length
                query_pos += length
            elif operation == 'D':
                for i in range(target_pos, target_pos + length):
                    freq[i][0]['D'] += 1
                    if i == 177637:
                        ovlp_pos.add((overlap.qname, None))

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
            if c1 < MIN_BASE:
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
    for tname, tovlps in overlaps.items():
        for overlap in tovlps:
            cigar = calculate_path(overlap, reads[tname], reads[overlap.qname])
            if cigar is not None:
                overlap.cigar = cigar

        #if tname == '34c1e790-6b60-4dc1-8160-0cbdb7ab85ff_1':
        #    print('Aln counts:', len(tovlps), count)


def trim_overlaps_cigar(
        cigar: List[Tuple[str, int]]) -> Tuple[int, int, int, int]:
    qstart_trim, tstart_trim = 0, 0
    while len(cigar) > 0 and (cigar[0][0] != '=' or cigar[0][1] < 5):
        op, length = cigar.pop(0)
        if op == '=' or op == 'X':
            qstart_trim += length
            tstart_trim += length
        elif op == 'I':
            qstart_trim += length
        elif op == 'D':
            tstart_trim += length
        else:
            ValueError('Invalid CIGAR operation.')

    qend_trim, tend_trim = 0, 0
    while len(cigar) > 0 and (cigar[-1][0] != '=' or cigar[-1][1] < 5):
        op, length = cigar.pop()
        if op == '=' or op == 'X':
            qend_trim += length
            tend_trim += length
        elif op == 'I':
            qend_trim += length
        elif op == 'D':
            tend_trim += length
        else:
            ValueError('Invalid CIGAR operation.')

    return qstart_trim, qend_trim, tstart_trim, tend_trim


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

    qstart_trim, qend_trim, tstart_trim, tend_trim = trim_overlaps_cigar(cigar)

    # Trim query
    if overlap.strand == '+':
        overlap.qstart += qstart_trim
        overlap.qend -= qend_trim
    else:
        overlap.qend -= qstart_trim
        overlap.qstart += qend_trim

    # Trim target
    overlap.tstart += tstart_trim
    overlap.tend -= tend_trim

    if len(cigar) == 0:
        return None

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


def median_n_overlaps(overlaps: Dict[str, List[Overlap]]) -> float:
    n_overlaps = [len(ovlps) for ovlps in overlaps.values()]
    return np.median(n_overlaps)


def main(args):
    # DEBUG_READ = 'eaf53782-10bf-45ab-a4f0-0052685ed651_1'

    true_overlaps = parse_paf(
        '/raid/scratch/stanojevicd/aisg/hg002/chr21/trimming/reads_to_ref/chr21_p_trimmed_to_ref.spaf'
    )

    overlaps = parse_paf(args.paf)
    # overlaps = {k: v for k, v in overlaps.items() if k == DEBUG_READ}
    overlaps = filter_primary_overlaps(overlaps)

    n_gt_removed, tn_removed, n_gt_small = 0, 0, 0
    valid_overlaps = {}
    for name, ovlps in overlaps.items():
        valid = []
        for o in ovlps:
            remove, tag = remove_overlap_ratio(o, true_overlaps)
            if tag == 0:
                n_gt_removed += 1
            if tag == 54:
                n_gt_small += 1
            if tag is None:
                tn_removed += 1
            if not remove:
                valid.append(o)

        if len(valid) > 0:
            valid_overlaps[name] = valid

    overlaps = valid_overlaps
    print('Gt overlaps removed because of ratio:', n_gt_removed, tn_removed,
          n_gt_small)

    extend_overlaps(overlaps)
    '''test_ovlps = set()
    with open('/home/stanojevicd/projects/EC_model/data/true_overlaps.txt',
              'r') as f:
        for line in f:
            line = line.strip()
            test_ovlps.add(line)

    count_in = 0
    for o in overlaps[DEBUG_READ]:
        if o.qname in test_ovlps:
            count_in += 1
    print('Number of true overlap', len(test_ovlps))
    print('Number of present true overlaps', count_in)'''
    valid_overlaps = {}
    n_valid, n_invalid = 0, 0
    counter = Counter()
    for name, ovlps in overlaps.items():
        valid = []
        # valid = [o for o in ovlps if not remove_overlap(o)]
        for o in ovlps:
            remove, tag = remove_overlap(o, true_overlaps)
            counter[tag] += 1
            if not remove:
                valid.append(o)

        if len(valid) > 0:
            valid_overlaps[name] = valid

        n_valid += len(valid)
        n_invalid += len(ovlps) - len(valid)
    print(f'Valid overlaps: {n_valid}, invalid overlaps: {n_invalid}')
    print(f'Counter {counter}')
    overlaps = valid_overlaps
    # sys.exit(-1)
    '''count_in = 0
    for o in overlaps[DEBUG_READ]:
        if o.qname in test_ovlps:
            count_in += 1
    print('Number of present true overlaps after filtering', count_in)'''

    #global median_ovlp_length
    # median_ovlp_length = median_overlap_length(overlaps)
    #print('Median overlap length:', median_ovlp_length)

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
            '''for r_id, record in reads.items():
                if r_id not in overlaps:
                    f_out.write(f'>{r_id} {record.description}\n')
                    f_out.write(f'{record.seq}\n')'''


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
    # N_OVERLAPS = int(1.2 * calculate_coverage(reads.values()) / 2)
    # N_OVERLAPS = 50

    t1 = time.time()
    main(args)
    t2 = time.time()
    print('Time elapsed: ', t2 - t1)
