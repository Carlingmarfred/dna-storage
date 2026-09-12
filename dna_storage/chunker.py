"""
Chunker module for DNA storage.
Splits data into fixed-size blocks with headers and reassembles them.
"""
import struct
import zlib
from typing import List


def chunk(data: bytes, block_sz: int = 256 * 1024) -> List[bytes]:
    """
    Splits data into fixed-size blocks with headers.
    
    Header format:
    - 4-byte block_id (uint32 BE)
    - 4-byte total_blocks (uint32 BE)
    - 4-byte payload_length (uint32 BE)
    - 4-byte CRC32 of the payload (uint32 BE)
    Total header size: 16 bytes.
    
    Args:
        data: The raw bytes to chunk.
        block_sz: The maximum size of the payload per block.
        
    Returns:
        A list of blocks, each containing a header and payload.
    """
    blocks = []
    data_len = len(data)
    
    if data_len == 0:
        return []
        
    total_blocks = (data_len + block_sz - 1) // block_sz
    
    for block_id in range(total_blocks):
        start_idx = block_id * block_sz
        end_idx = min(start_idx + block_sz, data_len)
        payload = data[start_idx:end_idx]
        
        payload_length = len(payload)
        crc32_val = zlib.crc32(payload) & 0xFFFFFFFF
        
        header = struct.pack('>IIII', block_id, total_blocks, payload_length, crc32_val)
        blocks.append(header + payload)
        
    return blocks


def validate_block(block: bytes) -> bool:
    """
    Checks CRC32 integrity of a single block.
    
    Args:
        block: The block bytes (header + payload).
        
    Returns:
        True if valid, False otherwise.
    """
    if len(block) < 16:
        return False
        
    header = block[:16]
    payload = block[16:]
    
    try:
        _, _, payload_length, expected_crc32 = struct.unpack('>IIII', header)
    except struct.error:
        return False
        
    if len(payload) != payload_length:
        return False
        
    actual_crc32 = zlib.crc32(payload) & 0xFFFFFFFF
    return actual_crc32 == expected_crc32


def reassemble(blocks: List[bytes]) -> bytes:
    """
    Re-assembles blocks in order, validates CRC32 of each block.
    
    Args:
        blocks: A list of block bytes.
        
    Returns:
        The reassembled raw bytes.
        
    Raises:
        ValueError: If a block is corrupted, missing, or invalid.
    """
    if not blocks:
        return b""
        
    parsed_blocks = []
    
    for block in blocks:
        if len(block) < 16:
            raise ValueError("Block too short to contain a header")
            
        header = block[:16]
        payload = block[16:]
        
        block_id, total_blocks, payload_length, expected_crc32 = struct.unpack('>IIII', header)
        
        if len(payload) != payload_length:
            raise ValueError(f"Block {block_id} has invalid payload length")
            
        actual_crc32 = zlib.crc32(payload) & 0xFFFFFFFF
        if actual_crc32 != expected_crc32:
            raise ValueError(f"Block {block_id} failed CRC32 validation")
            
        parsed_blocks.append((block_id, total_blocks, payload))
        
    # Sort blocks by block_id
    parsed_blocks.sort(key=lambda x: x[0])
    
    # Validation of sequence and total blocks
    if not parsed_blocks:
        return b""
        
    expected_total_blocks = parsed_blocks[0][1]
    if len(parsed_blocks) != expected_total_blocks:
        raise ValueError(f"Expected {expected_total_blocks} blocks, but got {len(parsed_blocks)}")
        
    reassembled_data = bytearray()
    for i, (block_id, total_blocks, payload) in enumerate(parsed_blocks):
        if block_id != i:
            raise ValueError(f"Missing block {i}")
        if total_blocks != expected_total_blocks:
            raise ValueError(f"Block {i} reports different total_blocks")
            
        reassembled_data.extend(payload)
        
    return bytes(reassembled_data)
