# DNA Storage Package
"""
🧬 dna_storage — AI-Optimized DNA Data Storage Pipeline
========================================================

Encode arbitrary digital payloads (audio, text, images) into synthetic DNA
oligos with AI-driven stability optimisation, Fountain + Reed-Solomon error
correction, and biological-constraint-aware nucleotide mapping.

Quick start
-----------
    >>> from dna_storage.binary_io import to_bytes, compress
    >>> from dna_storage.chunker import chunk, reassemble
    >>> from dna_storage.redundancy import HybridCodec
    >>> from dna_storage.dna_mapper import BalancedDNACodec
    >>> from dna_storage.ai_optimizer import DNAOptimizer
    >>> from dna_storage.primer_builder import PrimerBuilder
    >>> from dna_storage.synthesizer import OligoPool
    >>> from dna_storage.decoder import DNADecoder
    >>> from dna_storage.integrity import compute_sha256, verify_payload

Modules
-------
binary_io       Bytes ↔ file / string / JSON conversion + compression.
chunker         Fixed-size block splitting with CRC-32 headers.
redundancy      Fountain codes (LT) + Reed-Solomon hybrid error correction.
dna_mapper      Constraint-aware balanced 2-bit→base coding.
ai_optimizer    DNAStabilityTransformer + gradient-guided optimisation.
primer_builder  PCR primer / adapter attachment + barcoding.
synthesizer     Oligo pool → CSV/FASTA for synthesis providers.
decoder         Consensus builder + full decode pipeline.
integrity       CRC-32 / SHA-256 verification.
cli             Command-line interface (encode / decode / validate / info).
"""

__version__ = "1.0.0"
__author__ = "Genetic Data Layering Project"

from .binary_io import to_bytes, from_bytes, compress, decompress
from .chunker import chunk, reassemble, validate_block
from .integrity import compute_crc32, compute_sha256, verify_payload
