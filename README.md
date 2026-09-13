# DNA Storage

A pipeline for encoding digital data into synthetic DNA using error correction, fountain codes, sequence constraints, and AI-based optimization.

## Installation

```bash
git clone https://github.com/carlingmarfred/biological-nn
cd biological-nn
pip install -e .
```

## Usage

```bash
dna-storage encode --input data.zip --output out/
dna-storage decode --input out/oligos.fasta --output restored.zip
dna-storage simulate --input data.zip --error-rate 0.02
```

## Architecture

```text
Digital Data
     ↓
Error Correction
     ↓
Fountain Coding
     ↓
DNA Encoding
     ↓
AI Optimization
     ↓
Synthetic DNA
     ↓
Sequencing
     ↓
Decoding
     ↓
Original Data
```

## Core Components

* Reed-Solomon error correction
* Fountain coding
* Constraint-aware DNA encoding
* AI-based sequence optimization
* Consensus decoding
* Sequencing error simulation

## License

MIT
