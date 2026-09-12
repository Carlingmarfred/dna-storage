import pytest
import random
from dna_storage.chunker import chunk, reassemble, validate_block

def test_chunk_single_block():
    """Test chunking data smaller than the block size."""
    data = b'A' * 10
    blocks = chunk(data, block_sz=20)
    assert len(blocks) == 1
    assert reassemble(blocks) == data

def test_chunk_multiple_blocks():
    """Test chunking data larger than the block size."""
    data = b'A' * 25
    blocks = chunk(data, block_sz=10)
    assert len(blocks) == 3
    assert reassemble(blocks) == data

def test_chunk_exact_size():
    """Test chunking data exactly equal to the block size."""
    data = b'A' * 20
    blocks = chunk(data, block_sz=20)
    assert len(blocks) == 1
    assert reassemble(blocks) == data

def test_reassemble_in_order():
    """Test reassembling blocks in order."""
    data = b'Hello, world! This is a test.'
    blocks = chunk(data, block_sz=10)
    reassembled = reassemble(blocks)
    assert reassembled == data

def test_reassemble_out_of_order():
    """Test reassembling blocks that are shuffled."""
    data = b'Hello, world! This is a test.'
    blocks = chunk(data, block_sz=10)
    random.shuffle(blocks)
    reassembled = reassemble(blocks)
    assert reassembled == data

def test_validate_block_valid():
    """Test validating a valid block."""
    data = b'Valid data'
    blocks = chunk(data, block_sz=10)
    assert validate_block(blocks[0]) is True

def test_validate_block_corrupted():
    """Test validating a corrupted block."""
    data = b'Valid data'
    blocks = chunk(data, block_sz=10)
    block = bytearray(blocks[0])
    block[-1] ^= 0xFF
    assert validate_block(bytes(block)) is False

def test_roundtrip():
    """Test chunking and then reassembling returns original data."""
    data = bytes([random.randint(0, 255) for _ in range(100)])
    blocks = chunk(data, block_sz=15)
    assert reassemble(blocks) == data

def test_empty_data():
    """Test chunking empty data."""
    data = b''
    blocks = chunk(data, block_sz=10)
    assert len(blocks) == 0
    assert reassemble(blocks) == data
