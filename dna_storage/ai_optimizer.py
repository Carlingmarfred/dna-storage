"""
DNA stability optimizer module.

Uses a neural network to predict and improve the synthesizability of DNA oligos,
with a fallback to heuristic scoring if PyTorch is not available.
"""

import math
from typing import Dict, List, Tuple, Any

try:
    import torch
    import torch.nn as nn
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


def gc_content(seq: str) -> float:
    """Calculate the GC content of a DNA sequence."""
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for base in seq if base in 'GC')
    return gc / len(seq)


def max_homopolymer_run(seq: str) -> int:
    """Find the length of the longest homopolymer run in the sequence."""
    if not seq:
        return 0
    max_run = 1
    current_run = 1
    seq = seq.upper()
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            if current_run > max_run:
                max_run = current_run
        else:
            current_run = 1
    return max_run


def count_forbidden_motifs(seq: str) -> int:
    """Count the occurrences of forbidden restriction enzyme sites."""
    # EcoRI, BamHI, HindIII, NotI, XhoI
    motifs = ['GAATTC', 'GGATCC', 'AAGCTT', 'GCGGCCGC', 'CTCGAG']
    seq = seq.upper()
    count = 0
    for motif in motifs:
        # Simplistic non-overlapping count
        count += seq.count(motif)
    return count


def _reverse_complement(seq: str) -> str:
    """Get the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def predicted_hairpin_score(seq: str) -> float:
    """
    Calculate a simple self-complementarity score by checking if any 
    window of size 6-12 has its reverse complement appearing later in the sequence.
    Returns 1.0 if a large hairpin is found, down to 0.0 for none.
    """
    seq = seq.upper()
    max_hairpin_len = 0
    
    for window_size in range(12, 5, -1):
        if window_size > len(seq) // 2:
            continue
        found = False
        for i in range(len(seq) - window_size):
            window = seq[i:i + window_size]
            rc = _reverse_complement(window)
            if rc in seq[i + window_size:]:
                found = True
                break
        if found:
            max_hairpin_len = window_size
            break
            
    if max_hairpin_len >= 12:
        return 1.0
    elif max_hairpin_len <= 5:
        return 0.0
    else:
        return (max_hairpin_len - 5) / 7.0


def compute_features(seq: str) -> Dict[str, Any]:
    """Compute and return all sequence features as a dictionary."""
    return {
        'gc_content': gc_content(seq),
        'max_homopolymer_run': max_homopolymer_run(seq),
        'forbidden_motifs': count_forbidden_motifs(seq),
        'hairpin_score': predicted_hairpin_score(seq)
    }


class HeuristicStabilityScorer:
    """Scorer that evaluates DNA sequence stability using heuristic rules."""
    
    def score(self, seq: str) -> float:
        """
        Calculate a stability score in [0, 1].
        1.0 means perfectly stable/synthesizable.
        """
        features = compute_features(seq)
        
        # Penalties based on deviation and constraint violations
        gc_penalty = abs(features['gc_content'] - 0.5) * 2
        hp_penalty = features['max_homopolymer_run'] / max(1, len(seq))
        motif_penalty = features['forbidden_motifs'] * 0.1
        hairpin_penalty = features['hairpin_score'] * 0.15
        
        total_penalty = gc_penalty + hp_penalty + motif_penalty + hairpin_penalty
        
        score = 1.0 - total_penalty
        return max(0.0, min(1.0, score))


if TORCH_AVAILABLE:
    class PositionalEncoding(nn.Module):
        def __init__(self, d_model: int, max_len: int = 300):
            super().__init__()
            position = torch.arange(max_len).unsqueeze(1)
            div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
            pe = torch.zeros(max_len, 1, d_model)
            pe[:, 0, 0::2] = torch.sin(position * div_term)
            pe[:, 0, 1::2] = torch.cos(position * div_term)
            self.register_buffer('pe', pe)

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            x = x + self.pe[:x.size(0)]
            return x

    class DNAStabilityTransformer(nn.Module):
        """Transformer model to predict DNA sequence stability."""
        
        def __init__(self, d_model: int = 128, nhead: int = 8, dim_feedforward: int = 512, max_len: int = 300):
            super().__init__()
            self.input_projection = nn.Linear(4, d_model)
            self.pos_encoder = PositionalEncoding(d_model, max_len=max_len)
            
            encoder_layers = nn.TransformerEncoderLayer(
                d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward, batch_first=True
            )
            self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=6)
            
            self.output_head = nn.Sequential(
                nn.Linear(d_model, 1),
                nn.Sigmoid()
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            """
            Forward pass of the model.
            Args:
                x: one-hot encoded DNA sequence (batch, seq_len, 4)
            Returns:
                Stability scores, shape (batch, 1)
            """
            x = self.input_projection(x)
            
            # Positional encoding expects (seq_len, batch, d_model)
            x = x.transpose(0, 1)
            x = self.pos_encoder(x)
            x = x.transpose(0, 1)
            
            x = self.transformer_encoder(x)
            
            # Aggregate over sequence length
            x = x.mean(dim=1)
            
            # Output scalar [0, 1]
            scores = self.output_head(x)
            return scores


class DNAOptimizer:
    """Optimizer that improves DNA sequence stability iteratively."""
    
    def __init__(self, model=None, max_iterations: int = 10, use_heuristic: bool = True):
        self.model = model
        self.max_iterations = max_iterations
        self.use_heuristic = use_heuristic
        
        if self.model is None and self.use_heuristic:
            self.scorer = HeuristicStabilityScorer()
        elif self.model is not None and TORCH_AVAILABLE:
            self.scorer = self.model
        else:
            raise ValueError("No valid scoring method available.")

    def _score_sequence(self, seq: str) -> float:
        """Internal method to score a single sequence."""
        if hasattr(self.scorer, 'score'):
            return self.scorer.score(seq)
        elif TORCH_AVAILABLE and isinstance(self.scorer, nn.Module):
            self.scorer.eval()
            with torch.no_grad():
                mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
                tensor = torch.zeros((1, len(seq), 4), dtype=torch.float32)
                for i, base in enumerate(seq.upper()):
                    if base in mapping:
                        tensor[0, i, mapping[base]] = 1.0
                score = self.scorer(tensor)
                return score.item()
        return 0.0

    def optimize(self, seq: str) -> Tuple[str, float, Dict[str, Any]]:
        """
        Optimize a DNA sequence to improve its stability score.
        """
        current_seq = list(seq.upper())
        initial_score = self._score_sequence("".join(current_seq))
        current_score = initial_score
        iterations_used = 0
        substitutions_made = 0
        
        bases = ['A', 'C', 'G', 'T']
        
        for _ in range(self.max_iterations):
            best_seq = current_seq.copy()
            best_score = current_score
            improved = False
            
            for i in range(len(current_seq)):
                original_base = current_seq[i]
                for base in bases:
                    if base == original_base:
                        continue
                        
                    candidate_seq = current_seq.copy()
                    candidate_seq[i] = base
                    candidate_str = "".join(candidate_seq)
                    
                    score = self._score_sequence(candidate_str)
                    
                    if score > best_score:
                        best_score = score
                        best_seq = candidate_seq
                        improved = True
            
            if not improved:
                break
                
            current_seq = best_seq
            current_score = best_score
            iterations_used += 1
            substitutions_made += 1
            
        final_seq = "".join(current_seq)
        stats = {
            'iterations_used': iterations_used,
            'initial_score': initial_score,
            'final_score': current_score,
            'substitutions_made': substitutions_made
        }
        
        return final_seq, current_score, stats

    def optimize_batch(self, seqs: List[str]) -> List[Tuple[str, float, Dict[str, Any]]]:
        """Optimize a batch of sequences."""
        return [self.optimize(seq) for seq in seqs]


def save_model(model: Any, path: str) -> None:
    """Save the PyTorch model to a file."""
    if not TORCH_AVAILABLE:
        raise ImportError("PyTorch is not available. Cannot save model.")
    torch.save(model.state_dict(), path)


def load_model(path: str) -> Any:
    """Load the PyTorch model from a file."""
    if not TORCH_AVAILABLE:
        raise ImportError("PyTorch is not available. Cannot load model.")
    model = create_untrained_model()
    model.load_state_dict(torch.load(path, weights_only=True))
    return model


def create_untrained_model() -> Any:
    """Create a new, untrained instance of the DNAStabilityTransformer."""
    if not TORCH_AVAILABLE:
        raise ImportError("PyTorch is not available. Cannot create model.")
    return DNAStabilityTransformer()
