#!/usr/bin/env python3
import argparse


def fivePrimeFinder(pos: int, cigar: str, reverse: bool) -> int:
    """Return 5' coordinate accounting for strand and CIGAR."""
    cigar_ops = []
    num = ''
    for char in cigar:
        if char.isdigit():
            num += char
        else:
            cigar_ops.append((int(num), char))
            num = ''

    if not reverse:
        return pos 

    aligned_len = 0
    soft_clip_end = 0
    for i, (length, op) in enumerate(cigar_ops):
        if op in ('M', 'D', 'N'):
            aligned_len += length
        elif op == 'S' and i == len(cigar_ops) - 1:
            soft_clip_end = length
    return pos + aligned_len + soft_clip_end - 1


def isReverse(flag: int) -> bool:
    """
    Returns TRUE if the SAM FLAG indicates reverse strand
    """
    return (flag & 16) == 16

def getUMI(qname: str) -> str:
    """
    Extracts UMI from SAM QNAME field
    """
    return qname.split(":")[-1]

def dedup(SAM_in, SAM_out, UMI_file):
    seen_dict = {}
    prev_chrom = None
    dedup_count = 0
    invalid_umi_count = 0
    duplicate_count = 0

    with open(UMI_file) as f:
        umi_list = [line.strip() for line in f if line.strip()]

    with open(SAM_in) as in_sam, open(SAM_out, "w") as out_sam:
        for line in in_sam:
            if line.startswith("@"):
                out_sam.write(line)
                continue

            fields = line.split("\t")
            qname = fields[0]
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])
            cigar = fields[5]

            if chrom != prev_chrom:
                seen_dict = {}
                prev_chrom = chrom

            umi = getUMI(qname)
            if umi not in umi_list:
                invalid_umi_count += 1
                continue

            strand = isReverse(flag)
            five_pos = fivePrimeFinder(pos, cigar, strand)
            key = (chrom, five_pos, strand, umi)

            if key not in seen_dict:
                seen_dict[key] = True
                out_sam.write(line)
                dedup_count += 1
            else:
                duplicate_count += 1
    return dedup_count, duplicate_count, invalid_umi_count



def main():
    parser = argparse.ArgumentParser(description="Deduplication of single-end reads")
    parser.add_argument("-f", required=True, help="designates absolute file path to sorted sam file")
    parser.add_argument("-o", required=True, help="designates absolute file path to deduplicated sam file")
    parser.add_argument("-u", required=True, help="designates file containing unique molecular identifiers (UMIs)")
    args = parser.parse_args()

    unique, dupes, invalid = dedup(args.f, args.o, args.u)
    print("\nDedup complete.")
    print(f"Unique reads kept: {unique}")
    print(f"Duplicates removed: {dupes}")
    print(f"Invalid UMIs skipped: {invalid}")


if __name__ == "__main__":
    main()

