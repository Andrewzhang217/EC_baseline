from Bio import SeqIO
import edlib

from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path

from typing import *


@dataclass
class HAECSeqRecord:
    RC_BASE = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, name, id, description, seq):
        self.name = name
        self.id = id
        self.description = description
        self.seq = seq

    def reverse_complement(self,
                           start: Optional[int] = None,
                           end: Optional[int] = None) -> str:
        if start is None:
            start = 0
        if end is None:
            end = len(self.seq)

        result = ''.join(
            [self.RC_BASE[base] for base in reversed(self.seq[start:end])])
        return result


@dataclass
class Overlap:
    qname: str
    qstart: int
    qend: int
    tname: str
    tstart: int
    tend: int
    strand: str
    cigar: Optional[List[Tuple[str, int]]] = None


def get_reads(path: str) -> Dict[str, HAECSeqRecord]:
    file_type = Path(path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    read_records = {}
    for record in SeqIO.parse(path, file_type):
        read_records[record.id] = HAECSeqRecord(record.name, record.id,
                                                record.description,
                                                str(record.seq).upper())

    print("size: ", len(read_records))
    return read_records

def parse_paf(path: str, reads) -> Dict[str, List[Overlap]]:
    overlaps = defaultdict(list)
    counter = defaultdict(int)
    with open(path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            query_name = line[0]
            counter[query_name] += 1
    
    with open(path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')

            query_name = line[0]
            query_start = int(line[2])
            query_end = int(line[3])
            strand = line[4]
            target_name = line[5]
            target_start = int(line[7])
            target_end = int(line[8])

            overhang_threshold = 500
            query_len = len(reads[query_name].seq)
            target_len = len(reads[target_name].seq) 

            query_left = query_start - 0
            target_left = target_start - 0
            if counter[query_name] > 100 and query_left > overhang_threshold and target_left > overhang_threshold:
                if query_left < target_left:
                    query_head = reads[query_name].seq[0:100]
                    extend_start = target_start - query_left - 100
                    extend_end = target_start - query_left + 200
                    target_head = reads[target_name].seq[extend_start if extend_start > 0 else 0 : extend_end]
                    align = edlib.align(query_head, target_head, mode='HW')
                else :
                    target_head = reads[target_name].seq[0:100]
                    extend_start = query_start - target_left - 100
                    extend_end = query_start - target_left + 200
                    query_head = reads[query_name].seq[extend_start if extend_start > 0 else 0 : extend_end]
                    align = edlib.align(target_head, query_head, mode='HW')
                
                edit_distance = align['editDistance']
                if edit_distance > 15 :
                    continue

            query_right = query_len - query_end
            target_right = target_len - target_end
            if counter[query_name] > 100 and query_right > overhang_threshold and target_right > overhang_threshold:
                if query_right < target_right:
                    query_tail = reads[query_name].seq[query_len - 100:query_len]
                    extend_start = target_end + query_right - 200
                    extend_end = target_end + query_right + 100
                    target_tail = reads[target_name].seq[extend_start: extend_end if extend_end < target_len else target_len]
                    align = edlib.align(query_tail, target_tail, mode='HW')
                else :
                    target_tail = reads[target_name].seq[target_len - 100:target_len]
                    extend_start = query_end + target_right - 200
                    extend_end = query_end + target_right + 100
                    query_tail = reads[query_name].seq[extend_start: extend_end if extend_end < query_len else query_len]
                    align = edlib.align(target_tail, query_tail, mode='HW')
                edit_distance = align['editDistance']
                if edit_distance > 15 :
                    continue    

            # Save target
            overlap = Overlap(query_name, query_start, query_end, target_name,
                              target_start, target_end, strand)
            overlaps[target_name].append(overlap)

            # Save query as target
            overlap = Overlap(target_name, target_start, target_end, query_name,
                              query_start, query_end, strand)
            overlaps[query_name].append(overlap)

    # print('Before', len(overlaps['e012f204-6a49-4e82-884e-8138929a86c9_1']))
    return dict(overlaps)