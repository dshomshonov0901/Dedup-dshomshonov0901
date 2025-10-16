Problem: During RNA-seq library prep, PCR amplification increases the quantity of cDNA.
However, this also introduces PCR duplicates (repeated reads from the same molecule)
which will bias quanitification by overrepresenting transcripts.

Goal: Given a sorted SAM file of uniquely mapped single-end reads with known UMIs,
remove all PCR duplicates, keping only one copy per unique combination of:
Chromosome(RNAME) + Strand(FLAG) + 5' Start Position(POSITION) + UMI(QNAME)

Constraints: Millions of reads, account for soft clupping, ignore unknown UMIs, single-end only

Pseudocode:

1. Preprocessing:
Sort SAM file by chromosome and position using samtools:
samtools sort -O sam -O sorted_input.sam input.sam

2. Load 96 known UMIS into a set
   
3. Initialize:
   
-prev_chrom = NONE (keeps track of chromomsome of previous read processed)

-seen_dict = {} (stores every uniq combination of chrom, strand, start position, and umi. ex: seen_dict = {
   ('chr1', '+', 105, 'GAACAGGT'): True,
   ('chr1', '-', 206, 'TTGACCTA'): True
}

-open input_sam as read and output_sam as write

4. For each line in input_sam:
   
-If line starts with '@' -> write to output_sam (header)

-Else:

  a. Parse QNAME, FLAG, RNAME, POS, CIGAR
  
  b. Extract UMI from QNAME with getUMI function
      - If UMI is not in known set then skip read
      
  c. Determine strand: rev = TRUE if FLAG &16 == 16 else FALSE
  
  d. Calculate 5' position using fivePrimeFinder function
     - five_prime = fivePrimeFinder(POS, CIGAR, rev)
     
  e. key = (RNAME, rev, five_prime, UMI)
  
  f. if key not in seen_dict:
    - write line to output_sam
    - add key to seen dict
    
   Else:
      skip (this is duplicate!)
      
-If RNAME != prev_chrom:
  -Reset seen_dict
  -prev_chrom = RNAME

6. Close files

High Level Functions:

```python
def fivePrimeFinder(pos: int, cigar: str, reverse: bool) -> int:
    """
    Determines the 5′ start position for a read.
    Accounts for soft clipping, deletions, skipped regions, and strand.
    """
    for + strands: If there is soft clipping (S) at the start then subtract that length
    for - strands: Add up all reference consuming ops (M,N,D) + any soft clipping at the end
    Ignore insertions (I)

def isReverse(flag: int) -> bool:
    """
    Returns True if the SAM FLAG indicates a reverse strand.
    """
    return (flag & 16) == 16
    
def getUMI(qname: str) -> str:
    """
    Extracts the UMI from the SAM QNAME field.
    """
    return qname.split(":")[-1]

