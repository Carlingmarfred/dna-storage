"""
DNA mapping module for converting between binary data and DNA nucleotide sequences.
Enforces biological constraints like GC-balance, homopolymer limits, and avoiding restriction sites.
"""

import math
from typing import Dict, Any, List

BASES = ['A', 'C', 'G', 'T']
_BIT_TO_BASE = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
_BASE_TO_BIT = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}

def bytes_to_bits(data: bytes) -> str:
    """Convert bytes to a string of '0's and '1's."""
    return "".join(f"{b:08b}" for b in data)

def bits_to_bytes(bits: str) -> bytes:
    """Convert a string of '0's and '1's to bytes."""
    if len(bits) % 8 != 0:
        bits += '0' * (8 - len(bits) % 8)
    return bytes(int(bits[i:i+8], 2) for i in range(0, len(bits), 8))

def bits_to_dna_simple(data: bytes) -> str:
    """Simple 2-bit mapping from bytes to DNA string (00->A, 01->C, 10->G, 11->T)."""
    bits = bytes_to_bits(data)
    if len(bits) % 2 != 0:
        bits += '0'
    return "".join(_BIT_TO_BASE[bits[i:i+2]] for i in range(0, len(bits), 2))

def dna_to_bits_simple(seq: str) -> bytes:
    """Simple 2-bit mapping from DNA string to bytes."""
    bits = "".join(_BASE_TO_BIT[b] for b in seq)
    return bits_to_bytes(bits)

def complement(seq: str) -> str:
    """Return the DNA complement of the sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(comp[b] for b in seq)

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of the DNA sequence."""
    return complement(seq)[::-1]

def melting_temp_estimate(seq: str) -> float:
    """
    Estimate melting temperature (Tm) in Celsius.
    Uses Wallace rule for <=14 bases, simplified nearest neighbor/Marmur-Doty otherwise.
    """
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    at_count = seq.count('A') + seq.count('T')
    length = len(seq)
    
    if length == 0:
        return 0.0
    if length <= 14:
        return 2.0 * at_count + 4.0 * gc_count
    else:
        return 64.9 + 41.0 * (gc_count - 16.4) / length

class BalancedDNACodec:
    """
    Constraint-aware codec using a rotating cipher to prevent biological problems.
    Enforces GC balance, max homopolymer length, and avoids forbidden motifs.
    """
    
    FORBIDDEN_MOTIFS = ["GAATTC", "GGATCC", "AAGCTT"]  # EcoRI, BamHI, HindIII
    
    def __init__(self):
        pass
        
    @staticmethod
    def gc_content(seq: str) -> float:
        """Compute the fraction of G and C bases in the sequence."""
        if not seq:
            return 0.0
        return (seq.count('G') + seq.count('C')) / len(seq)
        
    @staticmethod
    def has_homopolymer(seq: str, max_run: int = 3) -> bool:
        """Check if the sequence contains runs of the same base longer than max_run."""
        if not seq:
            return False
        run = 1
        for i in range(1, len(seq)):
            if seq[i] == seq[i-1]:
                run += 1
                if run > max_run:
                    return True
            else:
                run = 1
        return False
        
    @staticmethod
    def has_forbidden_motif(seq: str) -> bool:
        """Check if the sequence contains any forbidden restriction enzyme sites."""
        for motif in BalancedDNACodec.FORBIDDEN_MOTIFS:
            if motif in seq:
                return True
        return False
        
    def validate(self, seq: str) -> Dict[str, Any]:
        """
        Validate the DNA sequence against all biological constraints.
        Returns a report dictionary.
        """
        gc = self.gc_content(seq)
        homo_viol = self.has_homopolymer(seq)
        motif_viol = self.has_forbidden_motif(seq)
        
        # Check windowed GC content (20 bases, 40-60%)
        window_gc_viol = False
        if len(seq) >= 20:
            for i in range(len(seq) - 19):
                window = seq[i:i+20]
                window_gc = self.gc_content(window)
                if window_gc < 0.4 or window_gc > 0.6:
                    window_gc_viol = True
                    break
        
        density = len(dna_to_bits_simple(seq)) * 8 / len(seq) if seq else 0.0
        
        return {
            "gc_content": gc,
            "window_gc_violation": window_gc_viol,
            "homopolymer_violation": homo_viol,
            "motif_violation": motif_viol,
            "length": len(seq),
            "density_bits_per_nt": density
        }
        
    def _get_alternative(self, base: str, index: int) -> str:
        """Get an alternative base mapping to avoid constraints."""
        # Simple rotation A->C->G->T->A
        alts = {'A': ['C', 'G', 'T'], 'C': ['G', 'T', 'A'], 'G': ['T', 'A', 'C'], 'T': ['A', 'C', 'G']}
        return alts[base][index % 3]

    def encode(self, data: bytes) -> str:
        """
        Convert bytes to constrained DNA using a constraint-aware alternating code.
        Bit 0 maps to 'A' (AT-rich) or 'C' (GC-rich).
        Bit 1 maps to 'T' (AT-rich) or 'G' (GC-rich).
        At each step, chooses the nucleotide that avoids repeating the previous base
        (preventing homopolymers) and keeps the overall GC content balanced near 50%.
        """
        bits = bytes_to_bits(data)
        seq: List[str] = []
        gc_count = 0
        prev = 'N'
        
        for bit in bits:
            if bit == '0':
                if prev == 'A':
                    choice = 'C'
                elif prev == 'C':
                    choice = 'A'
                else:
                    curr_gc = gc_count / len(seq) if seq else 0.5
                    choice = 'C' if curr_gc < 0.5 else 'A'
            else:
                if prev == 'G':
                    choice = 'T'
                elif prev == 'T':
                    choice = 'G'
                else:
                    curr_gc = gc_count / len(seq) if seq else 0.5
                    choice = 'G' if curr_gc < 0.5 else 'T'
                    
            if choice in ('G', 'C'):
                gc_count += 1
            seq.append(choice)
            prev = choice
            
        return "".join(seq)

    def decode(self, seq: str) -> bytes:
        """
        Reverse the constrained DNA mapping back to bytes.
        Bases 'A' and 'C' decode to bit 0.
        Bases 'G' and 'T' decode to bit 1.
        """
        bits: List[str] = []
        for b in seq:
            if b in ('A', 'C'):
                bits.append('0')
            elif b in ('G', 'T'):
                bits.append('1')
        bitstr = "".join(bits)
        return bits_to_bytes(bitstr)
