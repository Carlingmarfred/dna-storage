"""
Integrity verification module for checking payload integrity and computing hashes.
"""

import binascii
import hashlib
from typing import List, Dict, Any


def compute_crc32(data: bytes) -> int:
    """
    Compute the CRC32 checksum of the given data.

    Args:
        data: The binary data to checksum.

    Returns:
        The CRC32 integer value.
    """
    return binascii.crc32(data) & 0xFFFFFFFF


def compute_sha256(data: bytes) -> str:
    """
    Compute the SHA-256 hash of the given data.

    Args:
        data: The binary data to hash.

    Returns:
        The hex digest of the SHA-256 hash.
    """
    return hashlib.sha256(data).hexdigest()


def compute_md5(data: bytes) -> str:
    """
    Compute the MD5 hash of the given data.

    Args:
        data: The binary data to hash.

    Returns:
        The hex digest of the MD5 hash.
    """
    return hashlib.md5(data).hexdigest()


class IntegrityReport:
    """
    Report containing the results of payload integrity verification.
    """

    def __init__(self, original_hash: str, recovered_hash: str, block_results: List[Dict[str, Any]]):
        """
        Initialize the IntegrityReport.

        Args:
            original_hash: The expected hash of the original data.
            recovered_hash: The computed hash of the recovered data.
            block_results: List of dictionaries detailing verification for each block.
        """
        self.original_hash = original_hash
        self.recovered_hash = recovered_hash
        self.block_results = block_results

    @property
    def is_valid(self) -> bool:
        """Check if the recovered payload matches the original."""
        return self.original_hash == self.recovered_hash

    @property
    def corrupted_blocks(self) -> int:
        """Count the number of blocks that failed verification."""
        return sum(1 for res in self.block_results if not res.get("is_valid", True))

    @property
    def error_rate(self) -> float:
        """Calculate the ratio of corrupted blocks to total blocks."""
        if not self.block_results:
            return 0.0
        return self.corrupted_blocks / len(self.block_results)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the report to a dictionary representation.

        Returns:
            Dictionary containing report details.
        """
        return {
            "is_valid": self.is_valid,
            "original_hash": self.original_hash,
            "recovered_hash": self.recovered_hash,
            "total_blocks": len(self.block_results),
            "corrupted_blocks": self.corrupted_blocks,
            "error_rate": self.error_rate,
            "block_details": self.block_results
        }

    def __str__(self) -> str:
        """Get a human-readable summary of the report."""
        status = "PASSED" if self.is_valid else "FAILED"
        return (
            f"Integrity Verification: {status}\n"
            f"Original Hash (SHA-256): {self.original_hash}\n"
            f"Recovered Hash (SHA-256): {self.recovered_hash}\n"
            f"Blocks Corrupted: {self.corrupted_blocks} / {len(self.block_results)} "
            f"({self.error_rate * 100:.2f}%)"
        )


def verify_payload(original: bytes, recovered: bytes) -> IntegrityReport:
    """
    Verify the integrity of a recovered payload against the original.

    Args:
        original: The original binary data.
        recovered: The recovered binary data.

    Returns:
        An IntegrityReport summarizing the verification.
    """
    original_hash = compute_sha256(original)
    recovered_hash = compute_sha256(recovered)
    
    # In a full implementation, block_results would contain per-block CRC checks
    # Here we simulate with a general payload check as no blocks are provided
    is_valid = original_hash == recovered_hash
    block_results = [{"block_id": 0, "is_valid": is_valid}]

    return IntegrityReport(original_hash, recovered_hash, block_results)


def verify_block(block: bytes, expected_crc: int) -> bool:
    """
    Verify a single block using its expected CRC32 checksum.

    Args:
        block: The binary data of the block.
        expected_crc: The expected CRC32 integer value.

    Returns:
        True if the block's CRC32 matches the expected value, False otherwise.
    """
    return compute_crc32(block) == expected_crc
