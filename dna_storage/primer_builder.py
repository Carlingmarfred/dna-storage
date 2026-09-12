"""
Primer and adapter building for DNA storage sequences.
"""

import math
from typing import Dict, Any

ILLUMINA_FWD = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLUMINA_REV = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
UNIVERSAL_FWD = 'ATCGTACGATCGATCG'
UNIVERSAL_REV = 'CGATCGATCGTACGAT'

def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    mapping = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(mapping)[::-1]

def melting_temp(seq: str) -> float:
    """
    Calculates the basic melting temperature (Tm) of a sequence.
    Uses the Wallace rule for short sequences (<14 nt) and a
    simplified nearest-neighbor formulation for longer ones.
    """
    seq = seq.upper()
    g_count = seq.count('G')
    c_count = seq.count('C')
    a_count = seq.count('A')
    t_count = seq.count('T')
    
    if len(seq) < 14:
        return 2.0 * (a_count + t_count) + 4.0 * (g_count + c_count)
    else:
        return 64.9 + 41.0 * (g_count + c_count - 16.4) / len(seq)

def self_dimer_score(seq: str) -> float:
    """
    Calculates a basic self-dimer score by finding the length of the 
    longest reverse-complementary substring within the sequence.
    """
    max_len = 0
    seq = seq.upper()
    rc = reverse_complement(seq)
    # Simple longest common substring between seq and its reverse complement
    for i in range(len(seq)):
        for j in range(i + 1, len(seq) + 1):
            sub = seq[i:j]
            if sub in rc:
                max_len = max(max_len, len(sub))
    return float(max_len)

def generate_barcode(index_id: int, length: int = 8) -> str:
    """
    Converts an integer ID to a DNA barcode.
    Uses a simple base-4 encoding padded to the specified length.
    """
    bases = ['A', 'C', 'G', 'T']
    barcode = []
    val = index_id
    for _ in range(length):
        barcode.append(bases[val % 4])
        val //= 4
    return ''.join(barcode[::-1])

class PrimerBuilder:
    """
    Handles addition and validation of primers and indices for DNA storage.
    """
    def __init__(self, fwd_primer: str = UNIVERSAL_FWD, rev_primer: str = UNIVERSAL_REV, overlap_len: int = 20):
        self.fwd_primer = fwd_primer
        self.rev_primer = rev_primer
        self.overlap_len = overlap_len

    def add_primers(self, seq: str) -> str:
        """Prepends the forward primer and appends the reverse complement of the reverse primer."""
        return f"{self.fwd_primer}{seq}{reverse_complement(self.rev_primer)}"

    def strip_primers(self, seq: str) -> str:
        """Removes the forward and reverse primers from a sequenced read."""
        rc_rev = reverse_complement(self.rev_primer)
        start_idx = 0
        end_idx = len(seq)
        
        if seq.startswith(self.fwd_primer):
            start_idx = len(self.fwd_primer)
        
        if seq.endswith(rc_rev):
            end_idx = len(seq) - len(rc_rev)
            
        return seq[start_idx:end_idx]

    def add_index(self, seq: str, index_id: int, index_length: int = 8) -> str:
        """
        Adds a unique barcode index between the primer and payload for multiplexing.
        """
        barcode = generate_barcode(index_id, length=index_length)
        return f"{self.fwd_primer}{barcode}{seq}{reverse_complement(self.rev_primer)}"

    def validate_primers(self, seq: str) -> Dict[str, Any]:
        """
        Checks Tm compatibility, GC content, and self-dimer potential of the primers.
        """
        fwd_gc = (self.fwd_primer.count('G') + self.fwd_primer.count('C')) / len(self.fwd_primer) if self.fwd_primer else 0
        rev_gc = (self.rev_primer.count('G') + self.rev_primer.count('C')) / len(self.rev_primer) if self.rev_primer else 0
        
        return {
            'fwd_tm': melting_temp(self.fwd_primer),
            'rev_tm': melting_temp(self.rev_primer),
            'fwd_gc': fwd_gc,
            'rev_gc': rev_gc,
            'fwd_self_dimer': self_dimer_score(self.fwd_primer),
            'rev_self_dimer': self_dimer_score(self.rev_primer)
        }
