from Bio import SeqIO

from dataclasses import dataclass
from pathlib import Path

from typing import *


@dataclass
class HAECSeqRecord:
    RC_BASE = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, name: str, id: str, description: str, seq: str):
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

    return read_records