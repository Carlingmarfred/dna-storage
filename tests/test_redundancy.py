import pytest
from dna_storage.redundancy import (
    ReedSolomonCodec, 
    LTEncoder,
    LTDecoder,
    HybridCodec, 
    robust_soliton_distribution
)

def test_rs_roundtrip():
    """Test Reed-Solomon encode then decode returns original data."""
    codec = ReedSolomonCodec(nsym=10)
    data = b'Test RS encoding data.'
    encoded = codec.encode(data)
    decoded = codec.decode(encoded)
    assert decoded == data

def test_rs_error_correction():
    """Test Reed-Solomon recovers from introduced errors."""
    codec = ReedSolomonCodec(nsym=10)
    data = b'Test RS error correction.'
    encoded = bytearray(codec.encode(data))
    # Introduce a few errors
    encoded[0] = (encoded[0] + 1) % 256
    encoded[5] = (encoded[5] + 1) % 256
    decoded = codec.decode(bytes(encoded))
    assert decoded == data

def test_fountain_roundtrip():
    """Test LT encode and decode with enough packets recovers all blocks."""
    blocks = [b'BlockZero!', b'Block_One!', b'Block_Two!', b'BlockThree', b'Block_Four']
    encoder = LTEncoder(blocks, seed=42)
    packets = [encoder.generate_packet(i) for i in range(25)]
    
    decoder = LTDecoder(n_blocks=len(blocks), block_size=10)
    for pkt_id, indices, payload in packets:
        decoder.add_packet(pkt_id, indices, payload)
        if decoder.is_complete():
            break
            
    assert decoder.is_complete()
    assert decoder.get_blocks() == blocks

def test_fountain_insufficient_packets():
    """Test LT decoding fails gracefully with too few packets."""
    blocks = [b'A'*10, b'B'*10, b'C'*10, b'D'*10, b'E'*10]
    encoder = LTEncoder(blocks, seed=42)
    packets = [encoder.generate_packet(i) for i in range(2)]
    
    decoder = LTDecoder(n_blocks=len(blocks), block_size=10)
    for pkt_id, indices, payload in packets:
        decoder.add_packet(pkt_id, indices, payload)
        
    assert not decoder.is_complete()

def test_hybrid_roundtrip():
    """Test HybridCodec encode then decode returns original data."""
    codec = HybridCodec(rs_symbols=6, fountain_overhead=0.5, seed=123)
    data = b'Short test data for hybrid testing!' # exactly 35 bytes = 5 blocks of 7
    blocks = [data[i:i+7] for i in range(0, len(data), 7)]
    
    encoded = codec.encode(blocks)
    decoded = codec.decode(encoded, n_original_blocks=len(blocks), block_size=7)
    assert decoded == blocks

def test_soliton_distribution():
    """Test robust soliton distribution probabilities sum to ~1.0."""
    k = 100
    dist = robust_soliton_distribution(k, c=0.1, delta=0.5)
    total_prob = sum(dist)
    assert pytest.approx(total_prob, 0.01) == 1.0
