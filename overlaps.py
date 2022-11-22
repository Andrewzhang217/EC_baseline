from tqdm import tqdm

from dataclasses import dataclass
from collections import defaultdict
from math import ceil
import argparse
import sys

from sequences import *

from typing import *

COV_WINDOW_SIZE = 16
OL_THRESHOLD = 1
DOUBLE_OL_THRESHOLD = 2 * OL_THRESHOLD
EXTENSION_BP = 500


@dataclass
class Overlap:
    qname: str
    qlen: int
    qstart: int
    qend: int
    tname: str
    tlen: int
    tstart: int
    tend: int
    strand: str
    cigar: Optional[List[Tuple[str, int]]] = None


def parse_paf(path: str) -> Dict[str, List[Overlap]]:
    overlaps = defaultdict(list)

    with open(path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')

            query_name = line[0]
            query_length = int(line[1])
            query_start = int(line[2])
            query_end = int(line[3])
            strand = line[4]
            target_name = line[5]
            target_length = int(line[6])
            target_start = int(line[7])
            target_end = int(line[8])

            if query_name == target_name:
                continue

            # Save target
            overlap = Overlap(query_name, query_length, query_start, query_end,
                              target_name, target_length, target_start,
                              target_end, strand)
            overlaps[target_name].append(overlap)

            # Save query as target
            overlap = Overlap(overlap.tname, overlap.tlen, overlap.tstart,
                              overlap.tend, overlap.qname, overlap.qlen,
                              overlap.qstart, overlap.qend, overlap.strand)
            overlaps[query_name].append(overlap)

    return dict(overlaps)


def filter_primary_overlaps(
        overlaps: Dict[str, List[Overlap]]) -> Dict[str, List[Overlap]]:
    primary = {}
    for tname, ovlps in tqdm(overlaps.items(), 'Finding primary overlaps'):
        qname_to_ovlps = defaultdict(list)
        for o in ovlps:
            qname_to_ovlps[o.qname].append(o)

        prim_ovlps = [
            max(ovlps, key=lambda o: o.tend - o.tstart)
            for ovlps in qname_to_ovlps.values()
        ]
        primary[tname] = prim_ovlps

    return primary


def extend_overlaps(overlaps: Dict[str, List[Overlap]]) -> None:
    for ovlps in tqdm(overlaps.values(), 'Extending overlaps'):
        for o in ovlps:
            o.qstart = max(o.qstart - EXTENSION_BP, 0)
            o.qend = min(o.qend + EXTENSION_BP, o.qlen)
            o.tstart = max(o.tstart - EXTENSION_BP, 0)
            o.tend = min(o.tend + EXTENSION_BP, o.tlen)


def remove_overlap_ratio(overlap: Overlap, true_overlaps) -> bool:
    #t_ovlp_len = overlap.tend - overlap.tstart
    #q_ovlp_len = overlap.qend - overlap.qstart

    ratio = (overlap.tend - overlap.tstart) / (overlap.qend - overlap.qstart)

    if ratio > 1.111 or ratio < 0.9:
        '''try:
            for o in true_overlaps[overlap.tname]:
                if o.qname == overlap.qname:
                    print(overlap)
                    print(o)
                    print()
                    # print('Ratio filter\n')
                    if t_ovlp_len < 1000 or q_ovlp_len < 1000:
                        return True, 54
                    return True, 0
        except KeyError:
            pass'''

        return True, None  # Overlapping lengths differ significantly

    return False, -1


def remove_overlap(overlap: Overlap,
                   true_overlaps: Dict[str, List[Overlap]]) -> bool:
    if (overlap.qlen - (overlap.qend - overlap.qstart)) <= DOUBLE_OL_THRESHOLD:
        return False, 1  # Query covered completely

    if (overlap.tlen - (overlap.tend - overlap.tstart)) <= DOUBLE_OL_THRESHOLD:
        return False, 2  # Target covered completely

    # Fixing bug with wrongly filtering out prefix or suffix overlaps
    if overlap.strand == '-':
        qstart = overlap.qlen - overlap.qend
        qend = overlap.qlen - overlap.qstart
    else:
        qstart = overlap.qstart
        qend = overlap.qend

    if qstart > OL_THRESHOLD and overlap.tstart <= OL_THRESHOLD and (
            overlap.qlen - qend) <= OL_THRESHOLD:
        return False, 3  # Prefix overlap between query and target

    if overlap.tstart > OL_THRESHOLD and qstart <= OL_THRESHOLD and (
            overlap.tlen - overlap.tend) <= OL_THRESHOLD:
        return False, 4  # Suffix overlap between query and target

    try:
        for o in true_overlaps[overlap.tname]:
            if o.qname == overlap.qname:
                return True, 100
                '''print(overlap)
                print(o)
                print('Invalid overlap\n')'''
    except KeyError:
        pass

    return True, 5


def calc_pile(overlaps: List[Overlap], legnth: int) -> List[int]:
    cov = [0] * ceil(legnth / COV_WINDOW_SIZE)

    for o in overlaps:
        start_bin = o.tstart // COV_WINDOW_SIZE
        end_bin = ceil(o.tend / COV_WINDOW_SIZE)

        for i in range(start_bin, end_bin):
            cov[i] += 1

    return cov


def trim_and_split_reads(
        reads: Dict[str, HAECSeqRecord],
        overlaps: Dict[str, List[Overlap]]) -> List[HAECSeqRecord]:
    filtered_reads = []

    n_no_seq, n_multiple_regions, n_short = 0, 0, 0
    for r_id, record in tqdm(reads.items(), 'Trimming and splitting reads'):
        try:
            ovlps = overlaps[r_id]
        except KeyError:
            n_no_seq += 1
            continue

        pile = calc_pile(ovlps, len(record.seq))
        '''if r_id == '7bc08f84-931b-45c1-9d77-64b6b899bcf7_1':
            print(pile)
            sys.exit(-1)'''

        # Get start
        start = 0
        while start < len(pile) and pile[start] < 5:
            start += 1

        if start == len(pile):
            n_no_seq += 1
            continue

        end = len(pile) - 1  # Inclusive
        while end > start and pile[end] < 5:
            end -= 1
        end += 1  # Make it exclusive

        supported = True
        ranges = []
        for pos in range(start + 1, end):
            if supported and pile[pos] < 5:
                ranges.append((start, pos))
                supported = False
                continue

            if not supported and pile[pos] >= 5:
                start = pos
                supported = True
        ranges.append((start, end))

        if len(ranges) == 1:
            start, end = ranges[0]
            start, end = start * COV_WINDOW_SIZE, end * COV_WINDOW_SIZE  # Get exact positions

            if end - start >= 1000:
                record.seq = record.seq[start:end]
                filtered_reads.append(record)
            else:
                n_short += 1
        else:
            for i, (start, end) in enumerate(ranges):
                start, end = start * COV_WINDOW_SIZE, end * COV_WINDOW_SIZE

                if end - start >= 1000:
                    name = record.name + f'_{i}'
                    seq = record.seq[start:end]
                    filtered_reads.append(
                        HAECSeqRecord(name, name, record.description, seq))
                else:
                    n_short += 1

            n_multiple_regions += 1

    print('Number of unsupported reads:', n_no_seq)
    print('Number of reads with multiple supported regions:',
          n_multiple_regions)
    print('Number of removed segments (shorter than 500 bp):', n_short)

    return filtered_reads


def main(args: argparse.Namespace) -> None:
    reads = get_reads(args.reads)
    print('Number of reads:', len(reads))

    overlaps = parse_paf(args.paf)
    print('Number of overlaps:', sum([len(o) for o in overlaps.values()]))

    overlaps = filter_primary_overlaps(overlaps)
    reads = trim_and_split_reads(reads, overlaps)

    with open(args.output, 'w') as f:
        for read in reads:
            f.write(f'>{read.id} {read.description}\n')
            f.write(f'{read.seq}\n')
    print('Number of written reads:', len(reads))


def get_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--reads', type=str, required=True)
    parser.add_argument('-p', '--paf', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)

    return parser.parse_args()


if __name__ == '__main__':
    args = get_arguments()
    main(args)
