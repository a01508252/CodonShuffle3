# CodonShuffle3

### Description

This package is an updated version of the original [CodonShuffle](https://github.com/lauringlab/CodonShuffle/tree/master) tool, originally described by Jorge D M M, Mills R E, Lauring A S (2015). The updates ensure compatibility with Python3 and current computing environments, while preserving the core functionality: generating and analyzing permuted sequences of an open reading frame from a viral genome.

This updated package can be used to generate permuted sequences by shuffling bases in a way that preserves the protein coding sequence. The resulting sequences contain a large number of synonymous substitutions and may differ in various sequence-determined features (e.g., dinucleotide frequency or free energy of RNA folding). Additional scripts are then used to quantify these differences relative to the unpermuted sequence, and a least squares method is applied to identify permuted sequences that are most similar to the "wild type."

For more information or citation for the original CodonShuffle package, see:

* Jorge D M M, Mills R E, Lauring A S, CodonShuffle: a tool for generating and analyzing synonymously mutated sequences. [Virus Evolution](http://dx.doi.org/10.1093/ve/vev012), 2015, 1(1): vev012 

---

### Instructions

1. Download files and folders
2. Check software requirements
3. Install dependencies by running the install_dependencies.py script
4. Run the Python file (`CodonShuffle.py`)
5. Select sequence file in frame (`input_file.fas`)
6. Select a single permutation script to use (`n3`, `dn23`, `dn31`, `dn231`)
7. Select desired number of permuted sequences
8. Select the genomic feature to be evaluated

#### CodonShuffle.py

Script to shuffle nucleotides and evaluate genomic features

**Input Commands:**

- `-i` Input  
  *Fasta file (`input_file.fas`) of open reading frame beginning with ATG*

- `-s` Permutation script  
  *Choose one of: n3, dn23, dn31, dn231 to specify which permutation method to use*

- `-r` Number of replicates  
  *Number of permuted sequences that the program will generate, default is 1000*

- `-m` Genomic feature(s) to be used  
  *Specifies the genomic features for final least squares distance calculation. Defaults to all. To use a subset, list each feature: CAI, ENC, VFOLD, UFOLD, DN, CPB separated by spaces. If using RNA folding, specify algorithm (UNAfold or ViennaRNA) with this command.*

- `-g` Graphics  
  *Generates graphs of distributions of values for all genomic features*

- `--seed` Random seed  
  *Allows setting a random seed for reproducibility.*

- `-h` Help  
  *Displays help menu*

**Outputs:**

- Graphs of distribution of values for genomic features of permuted sequences. Each file will include the input sequence and permutation algorithm in its title, with suffixes such as:
    - `_dn.pdf` for dinucleotide frequency graphs
    - `_dnls.pdf` for overall dinucleotide bias least squares graphs
    - `fas.hamming.pdf` for Hamming distance graphs
    - `.fold.pdf` for RNA folding free energy graphs
    - `.out.enc.pdf` for effective number of codons (ENC) graphs
    - `.cai.pdf` for codon adaptation index (CAI) graphs
    - `.cpb.pdf` for codon pair bias (CPB) graphs

- A graph showing Hamming distance versus least squares distance: `fasfinal_graph.pdf`
- A table (`_final_table.txt`) with hamming distance, genomic feature values, and aggregate least squares distance. Headers include:
    - Sequence number (with input sequence as 0)
    - Distance.ls (overall least squares distance)
    - Nucleotide_difference
    - CPB (codon pair bias)
    - DN_least_square (aggregate dinucleotide least squares value)
    - Individual dinucleotide frequencies (e.g., DN..AA.)
    - VFOLD.mfe (minimum free energy from Vienna RNA)
    - ENC (effective number of codons)
    - CAI (codon adaptive index)
- A multi-sequence FASTA file (`.fas`) with all permuted sequences labeled (e.g., replicate_1)
- Additional intermediate output files during the CodonShuffle run:
    - `_least_square.txt`
    - `.blk`
    - `.cpb`
    - `.dn`
    - `.dnls`
    - `fas_distance_table.txt`
    - `.fasfold_table_mfe.txt`
    - `.out`
    - `new_table_final_graph`
    - `.cai`
    - `.fold`

**Example Usage:**

```bash
$ python CodonShuffle.py -i Poliovirus_1_Mahoney_P1.fas -s dn23 -r 100 -m CAI ENC CPB -g
