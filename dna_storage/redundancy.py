"""
redundancy.py

Provides a comprehensive redundancy module for DNA data storage featuring TWO layers
of error correction:
1. Reed-Solomon (RS) Codec for block-level error correction.
2. Fountain Code (LT Codes) for packet-loss resilience.
"""

import math
import random
from typing import List, Tuple, Dict, Set

try:
    import reedsolo
    from reedsolo import RSCodec
    REEDSOLO_AVAILABLE = True
except ImportError:
    REEDSOLO_AVAILABLE = False


class ReedSolomonCodec:
    """Wrapper for Reed-Solomon encoding and decoding using the reedsolo library."""
    
    def __init__(self, nsym: int = 32):
        """
        Initialize the Reed-Solomon codec.
        
        Args:
            nsym (int): Number of error correction symbols. Default 32 allows
                        correction of up to 16 symbol errors.
        """
        if not REEDSOLO_AVAILABLE:
            raise ImportError(
                "reedsolo library is not installed. "
                "Please install it (e.g. pip install reedsolo) to use ReedSolomonCodec."
            )
        self.nsym = nsym
        self.codec = RSCodec(nsym)

    def encode(self, data: bytes) -> bytes:
        """
        Apply Reed-Solomon encoding to the data.
        
        Args:
            data (bytes): The input data block.
            
        Returns:
            bytes: The encoded data block with parity symbols appended.
        """
        return bytes(self.codec.encode(data))

    def decode(self, data: bytes) -> bytes:
        """
        Decode and correct the Reed-Solomon encoded data.
        
        Args:
            data (bytes): The data block with parity symbols.
            
        Returns:
            bytes: The corrected original data block.
        """
        # decode returns a tuple: (decoded_data, decoded_message, err_pos)
        return bytes(self.codec.decode(data)[0])


def robust_soliton_distribution(k: int, c: float = 0.1, delta: float = 0.5) -> List[float]:
    """
    Generate the Robust Soliton distribution probabilities for LT codes.
    
    Args:
        k (int): Number of source blocks.
        c (float): Constant for the robust part of the distribution.
        delta (float): Probability of decoding failure.
        
    Returns:
        List[float]: A list of probabilities where index i corresponds to degree i+1.
    """
    if k <= 0:
        return []
    if k == 1:
        return [1.0]

    # Ideal soliton distribution
    rho = [0.0] * (k + 1)
    rho[1] = 1.0 / k
    for d in range(2, k + 1):
        rho[d] = 1.0 / (d * (d - 1))
        
    # Robust additions
    R = c * math.log(k / delta) * math.sqrt(k)
    tau = [0.0] * (k + 1)
    
    limit = int(round(k / R))
    if limit > k:
        limit = k
        
    for d in range(1, limit):
        tau[d] = R / (d * k)
    if limit <= k:
        tau[limit] = R * math.log(R / delta) / k
        
    # Combine and normalize
    mu = [0.0] * (k + 1)
    Z = 0.0
    for d in range(1, k + 1):
        mu[d] = rho[d] + tau[d]
        Z += mu[d]
        
    # Return probabilities for degrees 1 to k (0-indexed list)
    return [mu[d] / Z for d in range(1, k + 1)]


class LTEncoder:
    """Fountain Code Encoder using Luby Transform (LT) codes."""
    
    def __init__(self, blocks: List[bytes], seed: int = 42):
        """
        Initialize the LT Encoder.
        
        Args:
            blocks (List[bytes]): The input source blocks to encode.
            seed (int): Base random seed for reproducible packet generation.
        """
        self.blocks = blocks
        self.k = len(blocks)
        self.seed = seed
        self.distribution = robust_soliton_distribution(self.k) if self.k > 0 else []

    def generate_packet(self, packet_id: int) -> Tuple[int, List[int], bytes]:
        """
        Generate a single fountain packet.
        
        Args:
            packet_id (int): The identifier for the generated packet.
            
        Returns:
            Tuple[int, List[int], bytes]: A tuple containing the packet_id, 
                                          the list of source block indices XORed, 
                                          and the resulting XOR payload.
        """
        if self.k == 0:
            raise ValueError("No blocks to encode.")
            
        # Seed RNG with base seed + packet_id to ensure deterministic sampling
        local_rng = random.Random(self.seed + packet_id)
        
        # Sample the degree d based on the distribution
        r = local_rng.random()
        cumulative = 0.0
        d = self.k
        for i, p in enumerate(self.distribution):
            cumulative += p
            if r <= cumulative:
                d = i + 1
                break
                
        # Randomly choose d distinct blocks
        indices = local_rng.sample(range(self.k), d)
        
        # Compute the XOR payload of the chosen blocks
        payload = bytearray(self.blocks[indices[0]])
        for idx in indices[1:]:
            block = self.blocks[idx]
            for i in range(len(payload)):
                payload[i] ^= block[i]
                
        return (packet_id, indices, bytes(payload))

    def generate_packets(self, n_packets: int) -> List[Tuple[int, List[int], bytes]]:
        """
        Generate multiple fountain packets.
        
        Args:
            n_packets (int): Number of packets to generate.
            
        Returns:
            List[Tuple]: A list of generated packets.
        """
        return [self.generate_packet(i) for i in range(n_packets)]


class LTDecoder:
    """Fountain Code Decoder using belief propagation / peeling decoder."""
    
    def __init__(self, n_blocks: int, block_size: int):
        """
        Initialize the LT Decoder.
        
        Args:
            n_blocks (int): Total number of source blocks to recover.
            block_size (int): Size of each source block in bytes.
        """
        self.k = n_blocks
        self.block_size = block_size
        self.solved_blocks: Dict[int, bytearray] = {}
        
        # Maps packet_id -> {'indices': set of unsolved block indices, 'payload': current XORed bytes}
        self.unsolved_packets: Dict[int, dict] = {}
        
    def add_packet(self, packet_id: int, block_indices: List[int], payload: bytes) -> bool:
        """
        Incorporate a new received packet and run peeling decoding.
        
        Args:
            packet_id (int): The ID of the packet.
            block_indices (List[int]): Source block indices that were XORed.
            payload (bytes): The packet data payload.
            
        Returns:
            bool: True if all original blocks are completely recovered, False otherwise.
        """
        if self.is_complete():
            return True
            
        unsolved_indices: Set[int] = set()
        current_payload = bytearray(payload)
        
        # Pre-process packet: XOR out already solved blocks
        for idx in block_indices:
            if idx in self.solved_blocks:
                solved_block = self.solved_blocks[idx]
                for i in range(self.block_size):
                    current_payload[i] ^= solved_block[i]
            else:
                unsolved_indices.add(idx)
                
        # If no unsolved indices remain, the packet is redundant
        if len(unsolved_indices) == 0:
            return self.is_complete()
            
        # If exactly 1 unsolved index remains, we solve a new block!
        if len(unsolved_indices) == 1:
            idx = unsolved_indices.pop()
            self.solved_blocks[idx] = current_payload
            self._propagate()
        else:
            # Otherwise, store it for later
            self.unsolved_packets[packet_id] = {
                'indices': unsolved_indices,
                'payload': current_payload
            }
            
        return self.is_complete()
        
    def _propagate(self) -> None:
        """Run belief propagation across all unsolved packets."""
        progress = True
        while progress:
            progress = False
            
            # Iterate over a list of keys since we will modify the dictionary
            for pid in list(self.unsolved_packets.keys()):
                packet = self.unsolved_packets[pid]
                indices = packet['indices']
                payload = packet['payload']
                
                # XOR out newly solved blocks from this packet
                to_remove = set()
                for idx in indices:
                    if idx in self.solved_blocks:
                        solved_block = self.solved_blocks[idx]
                        for i in range(self.block_size):
                            payload[i] ^= solved_block[i]
                        to_remove.add(idx)
                        
                indices -= to_remove
                
                # Check if it has become solvable
                if len(indices) == 1:
                    idx = indices.pop()
                    if idx not in self.solved_blocks:
                        self.solved_blocks[idx] = bytearray(payload)
                        progress = True
                    del self.unsolved_packets[pid]
                elif len(indices) == 0:
                    del self.unsolved_packets[pid]

    def is_complete(self) -> bool:
        """Check if all source blocks have been decoded."""
        return len(self.solved_blocks) == self.k
        
    def get_blocks(self) -> List[bytes]:
        """
        Retrieve the decoded source blocks.
        
        Returns:
            List[bytes]: A list of the recovered original source blocks.
            
        Raises:
            ValueError: If the blocks are not completely decoded yet.
        """
        if not self.is_complete():
            raise ValueError(
                f"Decoding is not yet complete. Solved {len(self.solved_blocks)}/{self.k} blocks."
            )
        return [bytes(self.solved_blocks[i]) for i in range(self.k)]


class HybridCodec:
    """
    Hybrid Codec combining Reed-Solomon and LT Fountain Codes.
    Provides robust block-level corruption resilience via RS and 
    packet-loss resilience via Fountain Codes.
    """
    
    def __init__(self, rs_symbols: int = 32, fountain_overhead: float = 0.15, seed: int = 42):
        """
        Initialize the Hybrid Codec.
        
        Args:
            rs_symbols (int): Number of RS parity symbols per block.
            fountain_overhead (float): Extra fraction of LT packets to generate relative to block count.
            seed (int): Seed used by the fountain encoder.
        """
        self.rs_symbols = rs_symbols
        self.rs_codec = ReedSolomonCodec(nsym=rs_symbols)
        self.fountain_overhead = fountain_overhead
        self.seed = seed

    def encode(self, blocks: List[bytes]) -> List[Tuple[int, List[int], bytes]]:
        """
        Encode blocks by first applying RS encoding, then generating fountain packets.
        
        Args:
            blocks (List[bytes]): List of original data blocks.
            
        Returns:
            List[Tuple]: List of generated LT packets.
                         Note: Type is list of tuples, although prompt mentions list[bytes], 
                         this aligns structurally with the decoder requirement.
        """
        if not blocks:
            return []
            
        # 1. Reed-Solomon encode each block
        rs_blocks = [self.rs_codec.encode(b) for b in blocks]
        
        # 2. LT (Fountain) Encode
        k = len(rs_blocks)
        n_packets = int(math.ceil(k * (1.0 + self.fountain_overhead)))
        
        lt_encoder = LTEncoder(rs_blocks, seed=self.seed)
        packets = lt_encoder.generate_packets(n_packets)
        
        return packets

    def decode(self, packets: List[Tuple[int, List[int], bytes]], n_original_blocks: int, block_size: int) -> List[bytes]:
        """
        Decode packets by first resolving LT fountain code, then correcting with RS.
        
        Args:
            packets (List[Tuple]): List of received fountain packets.
            n_original_blocks (int): Number of blocks expected.
            block_size (int): Size of the ORIGINAL block before RS parity was added.
            
        Returns:
            List[bytes]: The fully decoded and corrected list of source blocks.
        """
        if not packets or n_original_blocks <= 0:
            return []
            
        # The LT decoder needs the size of the blocks it is reconstructing,
        # which include the RS parity symbols.
        rs_encoded_size = block_size + self.rs_symbols
        
        # 1. LT (Fountain) Decode
        lt_decoder = LTDecoder(n_blocks=n_original_blocks, block_size=rs_encoded_size)
        
        for pkt in packets:
            packet_id, indices, payload = pkt
            lt_decoder.add_packet(packet_id, indices, payload)
            if lt_decoder.is_complete():
                break
                
        if not lt_decoder.is_complete():
            raise ValueError(
                f"Fountain decoder failed to recover all blocks. "
                f"Solved {len(lt_decoder.solved_blocks)}/{n_original_blocks}."
            )
            
        rs_blocks = lt_decoder.get_blocks()
        
        # 2. Reed-Solomon Decode
        decoded_blocks = []
        for b in rs_blocks:
            decoded_blocks.append(self.rs_codec.decode(b))
            
        return decoded_blocks
