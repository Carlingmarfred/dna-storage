"""
Binary Input/Output utilities for DNA storage.
Handles conversions to/from bytes and basic compression.
"""
import base64
import json
import zlib
import struct
from pathlib import Path
from typing import Any, Union


def to_bytes(data: Union[str, Path, bytes, Any]) -> bytes:
    """
    Converts various data types into bytes.
    
    Args:
        data: The data to convert. Can be a string, a pathlib.Path (reads the file),
              raw bytes, or a JSON-serializable object.
              
    Returns:
        The byte representation of the data.
    """
    if isinstance(data, bytes):
        return data
    elif isinstance(data, Path):
        return data.read_bytes()
    elif isinstance(data, str):
        # Check if string happens to be a valid file path, if so read it
        try:
            path = Path(data)
            if path.is_file():
                return path.read_bytes()
        except OSError:
            pass
        return data.encode('utf-8')
    else:
        # Fallback to json serialization
        return json.dumps(data).encode('utf-8')


def from_bytes(blob: bytes, as_type: str) -> Any:
    """
    Converts bytes back to a specified type.
    
    Args:
        blob: The raw bytes.
        as_type: String specifying the return type ('str', 'json', 'base64', 'bytes').
        
    Returns:
        The converted data.
    """
    if as_type == 'bytes':
        return blob
    elif as_type == 'str':
        return blob.decode('utf-8')
    elif as_type == 'base64':
        return base64.b64encode(blob).decode('utf-8')
    elif as_type == 'json':
        return json.loads(blob.decode('utf-8'))
    else:
        raise ValueError(f"Unsupported type: {as_type}")


def compress(data: bytes) -> bytes:
    """
    Compresses data using zlib. 
    Adds a 2-byte header indicating compression status.
    Header: 0x00 0x01 (compressed) or 0x00 0x00 (uncompressed).
    
    Args:
        data: The raw bytes to compress.
        
    Returns:
        The compressed bytes with header.
    """
    try:
        compressed = zlib.compress(data)
        if len(compressed) < len(data):
            # 1 means compressed
            return struct.pack('>H', 1) + compressed
    except Exception:
        pass
    
    # 0 means uncompressed
    return struct.pack('>H', 0) + data


def decompress(data: bytes) -> bytes:
    """
    Decompresses data using zlib, checking the 2-byte header.
    
    Args:
        data: The compressed or uncompressed bytes with header.
        
    Returns:
        The original uncompressed bytes.
    """
    if len(data) < 2:
        return data
        
    header = struct.unpack('>H', data[:2])[0]
    payload = data[2:]
    
    if header == 1:
        return zlib.decompress(payload)
    else:
        return payload
