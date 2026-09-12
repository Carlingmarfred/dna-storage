import pytest
from dna_storage.dna_mapper import (
    bits_to_dna_simple, 
    dna_to_bits_simple, 
    BalancedDNACodec, 
    complement, 
    reverse_complement, 
    melting_temp_estimate
)

def test_simple_roundtrip():
    """Test bits_to_dna_simple and dna_to_bits_simple."""
    data = b'Hello'
    dna = bits_to_dna_simple(data)
    decoded_data = dna_to_bits_simple(dna)
    assert decoded_data == data
    
def test_balanced_codec_roundtrip():
    """Test encode and decode with BalancedDNACodec."""
    codec = BalancedDNACodec()
    data = b'Hello, world!'
    dna = codec.encode(data)
    decoded_data = codec.decode(dna)
    assert decoded_data == data

def test_gc_content_range():
    """Test that encoded sequences have a GC content between 0.2 and 0.8."""
    codec = BalancedDNACodec()
    data = b'\xff' * 50 + b'\x00' * 50 # Extreme data
    dna = codec.encode(data)
    gc_content = codec.gc_content(dna)
    assert 0.2 <= gc_content <= 0.8

def test_no_long_homopolymers():
    """Test that encoded sequences do not contain homopolymers > 3."""
    codec = BalancedDNACodec()
    data = b'\xff' * 50 + b'\x00' * 50
    dna = codec.encode(data)
    for base in 'ATCG':
        assert base * 4 not in dna

def test_complement():
    """Test DNA complement calculation."""
    assert complement('ATCG') == 'TAGC'

def test_reverse_complement():
    """Test DNA reverse complement calculation."""
    assert reverse_complement('ATCG') == 'CGAT'

def test_melting_temp():
    """Test melting temperature calculation returns reasonable values."""
    tm = melting_temp_estimate('ATCGATCGATCGATCG')
    assert 40 <= tm <= 80

def test_validate_report():
    """Test validate returns expected keys."""
    codec = BalancedDNACodec()
    report = codec.validate('ATCGATCG')
    assert 'gc_content' in report
    assert 'homopolymer_violation' in report
    assert 'motif_violation' in report
