"""
Synthesis order generator module for DNA storage.
"""

import csv
import statistics
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Any, Optional

@dataclass
class OligoPool:
    """
    Manages a pool of oligonucleotides for synthesis and formatting for providers.
    """
    oligos: List[str]
    metadata: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

    def to_csv(self, path: str, provider: str = 'generic') -> None:
        """
        Writes a CSV formatted for specific synthesis providers.
        """
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            
            if provider == 'generic':
                writer.writerow(['Name', 'Sequence', 'Length', 'GC_Content'])
                for i, seq in enumerate(self.oligos):
                    gc = (seq.count('G') + seq.count('C')) / len(seq) if seq else 0
                    writer.writerow([f'Oligo_{i}', seq, len(seq), f"{gc:.2f}"])
            
            elif provider == 'twist':
                writer.writerow(['Name', 'Sequence', 'Scale', 'Purification'])
                for i, seq in enumerate(self.oligos):
                    writer.writerow([f'Oligo_{i}', seq, '100pmol', 'Standard'])
            
            elif provider == 'idt':
                writer.writerow(['Name', 'Sequence', 'Scale', 'Purification'])
                for i, seq in enumerate(self.oligos):
                    writer.writerow([f'Oligo_{i}', seq, '25nmol', 'Standard'])
            
            else:
                raise ValueError(f"Unknown provider format: {provider}")

    def to_fasta(self, path: str) -> None:
        """
        Writes the oligo pool as a FASTA file.
        """
        with open(path, 'w') as f:
            for i, seq in enumerate(self.oligos):
                f.write(f">Oligo_{i}\n{seq}\n")

    def summary(self) -> Dict[str, Any]:
        """
        Returns statistics about the oligo pool.
        """
        if not self.oligos:
            return {'n_oligos': 0}
            
        lengths = [len(seq) for seq in self.oligos]
        gcs = [(seq.count('G') + seq.count('C')) / len(seq) for seq in self.oligos if seq]
        
        avg_len = sum(lengths) / len(lengths)
        return {
            'n_oligos': len(self.oligos),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'mean_length': avg_len,
            'mean_gc': statistics.mean(gcs) if gcs else 0.0,
            'estimated_cost': estimate_synthesis_cost(len(self.oligos), int(avg_len), provider='twist')
        }

    def validate(self) -> List[Dict[str, Any]]:
        """
        Validates each oligo for length limits, GC bounds, and homopolymers.
        """
        issues = []
        for i, seq in enumerate(self.oligos):
            seq_issues = []
            
            if len(seq) < 50 or len(seq) > 300:
                seq_issues.append(f"Length {len(seq)} outside 50-300 nt bounds")
                
            if seq:
                gc = (seq.count('G') + seq.count('C')) / len(seq)
                if gc < 0.2 or gc > 0.8:
                    seq_issues.append(f"GC content {gc:.2f} outside 0.2-0.8 bounds")
                    
            for base in ['A', 'C', 'G', 'T']:
                if base * 6 in seq:
                    seq_issues.append(f"Homopolymer of {base} found")
                    
            if seq_issues:
                issues.append({
                    'index': i,
                    'sequence': seq,
                    'issues': seq_issues
                })
                
        return issues

def estimate_synthesis_cost(n_oligos: int, avg_length: int, provider: str = 'twist') -> float:
    """
    Approximates the synthesis cost based on published pricing.
    """
    if provider == 'twist':
        # Hypothetical Twist cost: $0.09 per base for pool synthesis
        return n_oligos * avg_length * 0.09
    elif provider == 'idt':
        # Hypothetical IDT plate cost: $0.15 per base
        return n_oligos * avg_length * 0.15
    elif provider == 'generic':
        return n_oligos * avg_length * 0.10
    else:
        return 0.0
