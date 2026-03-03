# Deduper – PCR Duplicate Removal for SAM Files

## Overview

This project creates a reference-based PCR duplicate removal tool for **coordinate-sorted, single-end SAM files** using known Unique Molecular Identifiers (UMIs).

The script removes PCR duplicates and a read is considered a duplicate if it shares all of these:

- Chromosome  
- 5′ coordinate  
- Strand  
- UMI  

If all four match a previously seen read, the read is removed.

---

## How It Works

### 5′ Coordinate Calculation

Correct 5′ coordinate calculation is essential for accurate deduplication and needs to consider soft clipping.

- **Forward strand (FLAG 0)**  
  5′ = POS − left soft clipping

- **Reverse strand (FLAG 16)**  
  5′ = aligned end position + right soft clipping

CIGAR parsing accounts for reference-consuming operations (`M`, `D`, `N`, `=`, `X`) when computing the aligned length.

### Duplicate Tracking

Reads are tracked using a dictionary keyed by: (chromosome, five_prime_position, strand, UMI)

For efficiency, stored keys are cleared whenever the chromosome changes (input needs to be coordinate-sorted for this to work).

### Performance

- UMIs are stored in a **set**.
- The script processes the SAM file line-by-line without loading it entirely into memory.

---

## Requirements

- Python 3
- Coordinate-sorted SAM file
- Text file containing valid UMIs (one per line)

---

## Usage

Arguments
-f : Path to coordinate-sorted input SAM
-o : Path to deduplicated output SAM
-u : Path to UMI list file
After completion, the script reports:
Unique reads kept
Duplicates removed
Invalid UMIs skipped

Example usage:

```bash
python3 shomshonov_deduper.py \
-f input.sam \
-o output.sam \
-u UMI-list.txt
```

## Testing

Custom test files were created to validate:
Exact duplicate removal
Forward-strand soft clipping 
Reverse-strand soft clipping
Strand-specific deduplication
Different UMIs at the same coordinate
Invalid UMI filtering
Duplicate handling across different chromosomes
testfiles_readme.md explains the test files in more detail.
