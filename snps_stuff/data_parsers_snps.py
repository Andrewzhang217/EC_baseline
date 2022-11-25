from Bio import SeqIO

from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
from typing import *
from align_info import align_info


@dataclass
class HAECSeqRecord:
    RC_BASE = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, name, id, description, seq):
        self.name = name
        self.id = id
        self.description = description
        self.seq = seq
        self.info = None 
        self.hap = None

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
    #iden: float = None
    cigar: Optional[List[Tuple[str, int]]] = None

# delete infos after this function 
def get_reads_snps(path: str, infos: Dict[str, align_info]=None) -> Dict[str, HAECSeqRecord]:
    file_type = Path(path).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    read_records = {}
    if infos:
        for record in SeqIO.parse(path, file_type):
            record = HAECSeqRecord(record.name, record.id,
                        record.description,
                        str(record.seq).upper())
            
            record.info = infos[record.id]
            record.hap = int(record.id[-1])
            read_records[record.id] = record
    else:
        for record in SeqIO.parse(path, file_type):
            record = HAECSeqRecord(record.name, record.id,
                        record.description,
                        str(record.seq).upper())
            read_records[record.id] = record
    print("size: ", len(read_records))
    return read_records


def parse_paf(path: str) -> Dict[str, List[Overlap]]:
    overlaps = defaultdict(list)

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

           

            # remove self-overlap
            if query_name == target_name:
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
    print(len(overlaps))
    return dict(overlaps)