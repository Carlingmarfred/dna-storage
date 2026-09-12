"""
Command-line interface entry point for the dna_storage package.
"""

import argparse
import sys
import os

# Using relative imports to access sibling modules
# from .decoder import DNADecoder
# from .integrity import compute_sha256

def handle_encode(args):
    """Handle the 'encode' sub-command."""
    print(f"Encoding '{args.input}' to '{args.output}'")
    print(f"Format: {args.format}")
    print(f"Redundancy: {args.redundancy}")
    print(f"Block size: {args.block_size}")
    print(f"Compression enabled: {args.compress}")
    # Encoding logic would go here

def handle_decode(args):
    """Handle the 'decode' sub-command."""
    print(f"Decoding FASTQ '{args.input}' to '{args.output}'")
    if args.n_blocks:
        print(f"Expected blocks: {args.n_blocks}")
    if args.block_size:
        print(f"Block size: {args.block_size}")
    # Decoding logic would go here

def handle_validate(args):
    """Handle the 'validate' sub-command."""
    print(f"Validating oligo pool from '{args.input}'")
    # Validation logic would go here

def handle_info(args):
    """Handle the 'info' sub-command."""
    print(f"Analyzing file '{args.input}' for encoding stats...")
    # Info logic would go here

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="DNA Storage Pipeline CLI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available sub-commands")
    subparsers.required = True

    # Encode command
    encode_parser = subparsers.add_parser("encode", help="Encode a file into a DNA oligo pool")
    encode_parser.add_argument("--input", required=True, help="Input file path to encode")
    encode_parser.add_argument("--output", required=True, help="Output CSV path for the oligo pool")
    encode_parser.add_argument("--format", default="generic", choices=["generic", "twist", "idt"], help="Synthesis provider format")
    encode_parser.add_argument("--redundancy", type=float, default=0.15, help="Error correction redundancy (0.0 to 1.0)")
    encode_parser.add_argument("--block-size", type=int, default=262144, help="Block size in bytes")
    encode_parser.add_argument("--compress", action="store_true", help="Enable data compression before encoding")
    encode_parser.set_defaults(func=handle_encode)

    # Decode command
    decode_parser = subparsers.add_parser("decode", help="Decode sequencing reads into a file")
    decode_parser.add_argument("--input", required=True, help="Input FASTQ path containing sequencing reads")
    decode_parser.add_argument("--output", required=True, help="Output file path for the decoded payload")
    decode_parser.add_argument("--n-blocks", type=int, help="Expected number of blocks")
    decode_parser.add_argument("--block-size", type=int, help="Expected block size")
    decode_parser.set_defaults(func=handle_decode)

    # Validate command
    validate_parser = subparsers.add_parser("validate", help="Validate oligo pool constraints")
    validate_parser.add_argument("--input", required=True, help="Path to CSV/FASTA file containing oligos")
    validate_parser.set_defaults(func=handle_validate)

    # Info command
    info_parser = subparsers.add_parser("info", help="Show encoding statistics without encoding")
    info_parser.add_argument("--input", required=True, help="Path to the file to analyze")
    info_parser.set_defaults(func=handle_info)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
