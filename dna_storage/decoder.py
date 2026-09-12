"""
DNA decoding pipeline module.

This module provides components for assembling consensus sequences from sequencing reads
and decoding the consensus DNA sequences back into binary data.
"""

from typing import List, Tuple, Optional, Dict
from collections import defaultdict

def parse_fastq(path: str) -> List[Tuple[str, str, List[int]]]:
    """
    Parse a FASTQ file and extract sequences and quality scores.

    Args:
        path: Path to the FASTQ file.

    Returns:
        List of tuples containing (sequence_name, sequence, quality_scores).
        Quality scores are converted from ASCII phred (offset 33) to integers.
    """
    records = []
    with open(path, 'r', encoding='utf-8') as f:
        while True:
            name_line = f.readline()
            if not name_line:
                break
            name = name_line.strip()[1:]  # Remove '@'
            sequence = f.readline().strip()
            _ = f.readline()  # '+' line
            qual_line = f.readline().strip()
            
            # Phred+33 quality scores
            quality_scores = [ord(char) - 33 for char in qual_line]
            records.append((name, sequence, quality_scores))
            
    return records


class ConsensusBuilder:
    """
    Builds a consensus sequence from multiple sequencing reads using
    weighted majority voting based on quality scores.
    """

    def __init__(self, min_coverage: int = 3):
        """
        Initialize the ConsensusBuilder.

        Args:
            min_coverage: Minimum read depth required at a position to make a confident call.
        """
        self.min_coverage = min_coverage
        self.reads: List[str] = []
        self.quality_scores: List[Optional[List[int]]] = []

    def add_read(self, read: str, quality_scores: Optional[List[int]] = None) -> None:
        """
        Add a sequencing read to the builder.

        Args:
            read: The DNA sequence string.
            quality_scores: Optional list of Phred quality scores corresponding to the bases.
        """
        self.reads.append(read)
        self.quality_scores.append(quality_scores)

    def build_consensus(self) -> str:
        """
        Build the consensus sequence via per-position weighted majority vote.

        Returns:
            The consensus DNA sequence string.
        """
        if not self.reads:
            return ""

        max_len = max(len(r) for r in self.reads)
        consensus = []

        for i in range(max_len):
            counts = defaultdict(float)
            for read, q_scores in zip(self.reads, self.quality_scores):
                if i < len(read):
                    base = read[i]
                    weight = q_scores[i] if q_scores and i < len(q_scores) else 1.0
                    counts[base] += weight
            
            if counts:
                best_base = max(counts.items(), key=lambda x: x[1])[0]
                consensus.append(best_base)

        return "".join(consensus)

    def coverage_at(self, pos: int) -> int:
        """
        Get the read depth (coverage) at a specific position.

        Args:
            pos: The 0-indexed position.

        Returns:
            The number of reads covering the given position.
        """
        return sum(1 for r in self.reads if pos < len(r))

    def confidence_at(self, pos: int) -> float:
        """
        Calculate the confidence of the consensus base at a specific position.

        Args:
            pos: The 0-indexed position.

        Returns:
            Confidence value between 0.0 and 1.0, representing the weight
            fraction of the consensus base compared to all bases at that position.
        """
        counts = defaultdict(float)
        total_weight = 0.0

        for read, q_scores in zip(self.reads, self.quality_scores):
            if pos < len(read):
                base = read[pos]
                weight = q_scores[pos] if q_scores and pos < len(q_scores) else 1.0
                counts[base] += weight
                total_weight += weight

        if total_weight == 0.0:
            return 0.0

        max_weight = max(counts.values()) if counts else 0.0
        return max_weight / total_weight


class DNADecoder:
    """
    Decodes sequenced DNA reads back into the original digital payload.
    """

    def __init__(self, primer_builder=None, balanced_codec=None, hybrid_codec=None):
        """
        Initialize the DNADecoder.

        Args:
            primer_builder: Component to handle primer stripping and indexing.
            balanced_codec: Component to convert DNA sequences to bits.
            hybrid_codec: Component to handle Fountain/Reed-Solomon decoding.
        """
        self.primer_builder = primer_builder
        self.balanced_codec = balanced_codec
        self.hybrid_codec = hybrid_codec

    def decode_reads(self, reads: List[str], quality_scores: Optional[List[List[int]]] = None) -> bytes:
        """
        Decode a set of raw sequencing reads into the original binary payload.

        Args:
            reads: List of DNA sequence strings.
            quality_scores: Optional list of quality score lists for each read.

        Returns:
            The decoded binary payload.
        """
        # Step 1: Strip primers and extract indices (mock logic if components missing)
        # Step 2: Group by barcode index
        groups: Dict[int, ConsensusBuilder] = defaultdict(ConsensusBuilder)
        
        for idx, read in enumerate(reads):
            q_scores = quality_scores[idx] if quality_scores else None
            
            if self.primer_builder:
                stripped, block_idx = self.primer_builder.strip_primers(read)
            else:
                # Mock fallback
                stripped, block_idx = read, 0

            groups[block_idx].add_read(stripped, q_scores)

        # Step 3: Build consensus per group
        consensus_sequences = {}
        for block_idx, builder in groups.items():
            consensus_sequences[block_idx] = builder.build_consensus()

        # Step 4: Decode DNA -> bits
        encoded_blocks = []
        for block_idx in sorted(consensus_sequences.keys()):
            dna_seq = consensus_sequences[block_idx]
            if self.balanced_codec:
                bits = self.balanced_codec.decode(dna_seq)
            else:
                bits = dna_seq.encode('utf-8')  # Mock
            encoded_blocks.append(bits)

        # Step 5: Apply hybrid codec (Fountain + RS)
        if self.hybrid_codec:
            original_blocks = self.hybrid_codec.decode(encoded_blocks)
        else:
            original_blocks = encoded_blocks

        # Step 6: Reassemble
        payload = b"".join(original_blocks)
        return payload

    def decode_fastq(self, fastq_path: str) -> bytes:
        """
        Parse a FASTQ file and decode its contents into the original payload.

        Args:
            fastq_path: Path to the FASTQ file.

        Returns:
            The decoded binary payload.
        """
        records = parse_fastq(fastq_path)
        reads = [r[1] for r in records]
        quality_scores = [r[2] for r in records]
        return self.decode_reads(reads, quality_scores)
