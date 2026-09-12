# 🧬 DNA Storage — AI-Optimized Digital Data Encoding in Synthetic DNA

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/carlingmarfred/biological-nn)

An advanced pipeline for encoding digital information into synthetic DNA. This project leverages fountain codes, Reed-Solomon error correction, constraint-aware balanced DNA coding, and optional neural network optimization to ensure stable, high-density, and reliable data storage in synthetic biological sequences.

---

## 🏗️ Architecture

```text
[ Digital File ] 
       ↓ 
[ Byte Chunker ] ------(Splits data into manageable chunks)
       ↓
[ Error Correction ] --(Adds Reed-Solomon parity bytes)
       ↓
[ Fountain Coder ] ----(Generates practically infinite droplets)
       ↓
[ DNA Encoder ] -------(Maps bytes to A,C,G,T adhering to GC% & homopolymer limits)
       ↓
[ AI Optimizer ] ------(Neural Network scoring & filtering for sequence stability)
       ↓
[ Payload Assembler ] -(Attaches addressing headers, primers, and adapters)
       ↓
[ Synthesized DNA Pools ]
       ↓
[ Sequencer ]
       ↓
[ Consensus Decoder ]
       ↓
[ Original File ]
```

---

## 🚀 Installation

Ensure you have Python 3.10+ installed.

### Basic Installation (Core Pipeline)
```bash
git clone https://github.com/yourusername/dna-storage.git
cd dna-storage
pip install -e .
```

### Full Installation (with AI optimization and Dev tools)
```bash
pip install -e .[all]
```

---

## ⚡ Quick Start

### Python API

```python
from dna_storage.pipeline import run_pipeline, decode_pipeline
import os

# 1. Encode a file to DNA
input_file = "my_data.txt"
with open(input_file, "w") as f:
    f.write("Hello, DNA world!")

encoded_sequences = run_pipeline(
    input_filepath=input_file,
    output_dir="./output",
    redundancy=0.2,
    use_ai=True
)

print(f"Generated {len(encoded_sequences)} oligos.")

# 2. Decode DNA back to file
decoded_file = decode_pipeline(
    input_fasta="./output/encoded_oligos.fasta",
    output_filepath="./decoded_data.txt"
)

with open(decoded_file, "r") as f:
    print(f"Decoded: {f.read()}")
```

---

## 💻 CLI Usage

The package provides a convenient command-line interface with 4 sub-commands:

- **`encode`**: Convert a digital file into a FASTA file containing synthetic DNA oligos.
- **`decode`**: Decode a FASTA or FASTQ file back into the original digital file.
- **`simulate`**: Run a full encode->error simulation->decode loop to test pipeline robustness.
- **`optimize`**: Train or tune the AI stability model on your hardware.

**Examples:**
```bash
# Encode with 30% redundancy
dna-storage encode --input data.zip --output out_dir --redundancy 0.3

# Decode sequences
dna-storage decode --input out_dir/oligos.fasta --output restored_data.zip

# Simulate the full pipeline with a 2% dropout rate
dna-storage simulate --input test.txt --error-rate 0.02
```

---

## 📚 Module Reference

| Module | Description |
|--------|-------------|
| `cli.py` | Command-line interface entry point. |
| `pipeline.py` | High-level orchestration of the encode and decode processes. |
| `io_utils.py` | File reading, writing, and FASTA/FASTQ handling. |
| `chunking.py` | Splitting payloads and managing addressing metadata. |
| `error_correction.py` | Reed-Solomon encoding and decoding algorithms. |
| `fountain.py` | Droplet generation and robust set covering for data recovery. |
| `dna_codec.py` | Constraint-aware mapping between binary data and ACGT nucleotides. |
| `ai_optimizer.py` | Neural network models to predict and optimize sequence stability. |
| `assembly.py` | Attachment and trimming of primers, adapters, and payload formatting. |
| `consensus.py` | Algorithms for clustering reads and generating consensus sequences. |
| `simulation.py` | Synthesizer and sequencer error simulation (mutations, dropouts). |

---

## ⚙️ How It Works

1. **Binary Conversion & Chunking**: The input file is converted to a binary stream and chunked into fixed-size blocks. Addressing metadata is attached to each chunk so its original position is known.
2. **Fountain Codes + Reed-Solomon Redundancy**: Reed-Solomon parity bytes are added to intra-chunk data to fix point mutations. A fountain coding scheme generates an unlimited supply of 'droplets' to ensure robustness against sequence dropout.
3. **Constraint-Aware Balanced DNA Coding**: Bytes are mapped to base pairs (A,C,G,T) using specialized dictionaries or state machines that strictly maintain 40-60% GC content and prevent homopolymer runs (>4 identical bases).
4. **AI Stability Optimization**: An optional neural network scores generated sequences for secondary structure formation and biochemical stability. Poor candidates are regenerated or altered.
5. **Primer/Adapter Attachment**: Constant region sequences are prepended and appended for PCR amplification and sequencer compatibility.
6. **Sequencing & Consensus Decoding**: After sequencing, reads are clustered and aligned. A consensus algorithm extracts the most likely original sequence, and the decoding layers (RS + Fountain) reconstruct the original file.

---

## 📊 Performance Metrics

| Metric | Target Specification |
|--------|----------------------|
| **Logical Density** | ~1.6 bits / nucleotide |
| **Recovery Rate** | 100% recovery with up to 15% missing sequences |
| **Encoding Speed** | >10 MB / minute (CPU) |
| **Decoding Speed** | >5 MB / minute (CPU) |
| **AI Optimization** | ~1M bases scored / second (GPU) |

*(Performance varies based on hardware and dataset size)*

---

## 🤝 Contributing

Contributions are welcome! Please open an issue to discuss proposed changes before submitting a Pull Request.
Ensure you run tests and maintain type hint coverage:
```bash
pytest tests/
mypy dna_storage/
```

---

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
