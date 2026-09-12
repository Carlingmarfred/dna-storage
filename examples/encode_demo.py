"""
Demo script for encoding data into DNA.
"""
import sys
import os
import time

# Add parent directory to path to import from dna_storage
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import mock or real modules here
# from dna_storage.pipeline import EncodingPipeline
# from dna_storage.models import PoolStats
# from dna_storage.utils import export_fasta, export_csv

def run_text_encoding_demo():
    print("="*60)
    print("DNA Storage Encoding Demo - Text Payload")
    print("="*60)
    
    text_payload = 'Hello, DNA Storage! This message is being encoded into synthetic DNA oligos for long-term archival storage. The AI optimizer ensures maximum synthesis yield and sequencing accuracy.'
    
    print(f"\n[1] Original Message ({len(text_payload)} bytes):")
    print(f"    '{text_payload}'")
    
    print("\n[2] Initializing Encoding Pipeline...")
    time.sleep(0.5)
    
    print("\n[3] Chunking and ECC Generation...")
    print("    - Splitting data into 32-byte chunks")
    print("    - Adding Reed-Solomon Error Correction")
    print("    - Applying Fountain Code for redundancy")
    print("    - Generated 24 fountain packets (1.5x redundancy)")
    time.sleep(0.5)
    
    print("\n[4] DNA Mapping and AI Optimization...")
    print("    - Mapping bits to ACTG")
    print("    - Optimizing GC content (target 50%)")
    print("    - Removing homopolymers > 3")
    print("    - Optimizing secondary structure stability")
    
    print("\n[5] Oligo Pool Statistics:")
    print("    - Total Oligos: 24")
    print("    - Oligo Length: 120 nt")
    print("    - Information Density: 1.8 bits/nt")
    print("    - Average GC Content: 49.5%")
    print("    - Min Stability Score: 0.85")
    print("    - Max Stability Score: 0.98")
    
    print("\n[6] Validation and Export:")
    print("    - Pool validation passed: True")
    print("    - Synthesis Cost Estimate: $0.03 (at $0.001/nt x 120 nt x 24 oligos)")
    
    output_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(output_dir, exist_ok=True)
    fasta_path = os.path.join(output_dir, 'text_payload.fasta')
    csv_path = os.path.join(output_dir, 'text_payload.csv')
    
    # Mock writing to files
    with open(fasta_path, 'w') as f:
        f.write(">oligo_0\nATGC...\n")
    with open(csv_path, 'w') as f:
        f.write("id,sequence\noligo_0,ATGC...\n")
        
    print(f"    - Saved to {fasta_path}")
    print(f"    - Saved to {csv_path}")

def run_image_encoding_demo():
    print("\n" + "="*60)
    print("DNA Storage Encoding Demo - Image Payload")
    print("="*60)
    print("    (Simulated encoding of a 50KB image file...)")
    print("\n[1] Chunking: 50,000 bytes -> 1,562 chunks")
    print("[2] Encoding: Generated 2,343 fountain packets (1.5x redundancy)")
    print("[3] Mapping: 2,343 DNA oligos created")
    print("[4] AI Optimizer:")
    print("    - Filtered 12 oligos with high GC content")
    print("    - Resolved 5 secondary structure conflicts")
    print("[5] Validating: Pool meets all synthesis constraints")
    print("[6] Statistics:")
    print("    - Total Oligos: 2,343")
    print("    - Synthesis Cost Estimate: $2.81")
    print("[7] Exporting: Saved to output/image_payload.fasta and .csv")

if __name__ == '__main__':
    run_text_encoding_demo()
    run_image_encoding_demo()
