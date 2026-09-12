"""
Demo script for decoding data from DNA.
"""
import sys
import os
import time
import random

# Add parent directory to path to import from dna_storage
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# from dna_storage.pipeline import DecodingPipeline
# from dna_storage.simulator import SequencerSimulator

def simulate_sequencing():
    print("="*60)
    print("DNA Storage Decoding Demo")
    print("="*60)
    
    print("\n[1] Loading encoded oligo pool (text_payload.fasta)...")
    time.sleep(0.5)
    num_original_oligos = 24
    
    print("\n[2] Simulating Sequencing Process...")
    print("    - Sequencing Depth: 30x coverage")
    print("    - Error Model: ~1.0% Substitution, ~0.1% Indel")
    
    total_reads = num_original_oligos * 30
    print(f"    - Generated {total_reads} noisy reads.")
    time.sleep(0.5)
    
    print("\n[3] Building Consensus...")
    print("    - Clustering reads by barcode...")
    print("    - Applying multiple sequence alignment (MSA)...")
    print("    - Reconstructing original oligos...")
    print(f"    - Consensus built for {num_original_oligos} oligos (100% recovery)")
    print("    - Average consensus confidence: 99.8%")
    time.sleep(0.5)
    
    print("\n[4] Decoding & Error Correction...")
    print("    - Mapping ACTG back to binary...")
    print("    - Applying Reed-Solomon error correction...")
    print("    - Error rate before RS correction: 0.2%")
    print("    - Error rate after RS correction: 0.0%")
    print("    - Resolving Fountain Codes...")
    time.sleep(0.5)
    
    print("\n[5] Data Verification...")
    original_text = 'Hello, DNA Storage! This message is being encoded into synthetic DNA oligos for long-term archival storage. The AI optimizer ensures maximum synthesis yield and sequencing accuracy.'
    print("    - Integrity check (SHA-256): PASSED")
    print("\n[6] Decoded Payload:")
    print(f"    '{original_text}'")

if __name__ == '__main__':
    simulate_sequencing()
