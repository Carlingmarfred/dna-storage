import pytest
import random
import hashlib
from dna_storage.chunker import chunk, reassemble
from dna_storage.redundancy import ReedSolomonCodec
from dna_storage.dna_mapper import BalancedDNACodec
from dna_storage.binary_io import compress, decompress

def encode_pipeline(data: bytes) -> list[str]:
    blocks = chunk(data, block_sz=16)
    rs_codec = ReedSolomonCodec(nsym=4)
    dna_codec = BalancedDNACodec()
    
    oligos = []
    for block in blocks:
        protected = rs_codec.encode(block)
        oligos.append(dna_codec.encode(protected))
    return oligos

def decode_pipeline(oligos: list[str]) -> bytes:
    dna_codec = BalancedDNACodec()
    rs_codec = ReedSolomonCodec(nsym=4)
    
    recovered_blocks = []
    for oligo in oligos:
        payload_bytes = dna_codec.decode(oligo)
        recovered_block = rs_codec.decode(payload_bytes)
        recovered_blocks.append(recovered_block)
        
    return reassemble(recovered_blocks)

def test_text_roundtrip():
    """Test full pipeline roundtrip for text data."""
    text = "This is a test of the DNA storage pipeline. It should encode and decode perfectly."
    data = text.encode('utf-8')
    oligos = encode_pipeline(data)
    decoded_data = decode_pipeline(oligos)
    assert decoded_data == data

def test_binary_roundtrip():
    """Test full pipeline roundtrip for binary data."""
    data = bytes([random.randint(0, 255) for _ in range(128)])
    oligos = encode_pipeline(data)
    decoded_data = decode_pipeline(oligos)
    assert decoded_data == data

def test_compressed_roundtrip():
    """Test full pipeline roundtrip with compression."""
    data = b"A" * 100 + b"B" * 100
    compressed = compress(data)
    oligos = encode_pipeline(compressed)
    decoded_data = decode_pipeline(oligos)
    decompressed = decompress(decoded_data)
    assert decompressed == data

def test_integrity_verification():
    """Verify SHA-256 matches before and after roundtrip."""
    data = bytes([random.randint(0, 255) for _ in range(256)])
    original_hash = hashlib.sha256(data).hexdigest()
    
    oligos = encode_pipeline(data)
    decoded_data = decode_pipeline(oligos)
    
    recovered_hash = hashlib.sha256(decoded_data).hexdigest()
    assert original_hash == recovered_hash
