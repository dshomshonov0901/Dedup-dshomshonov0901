Test cases covered in input_sam

1) Exact duplicate removal (forward strand, same chr/pos/cigar/UMI):
   - Reads :1139 and :1140 are identical keys so only :1139 should remain.

2) Forward-strand soft-clip affects 5′ coordinate (dedup via adjusted 5′):
   - :1141 is chr1 POS=95 CIGAR=50M. 5′ = 95
   - :1142 is chr1 POS=100 CIGAR=5S45M. 5′ = 100 - 5 = 95
   Same chr/strand/UMI and same adjusted 5′ so :1142 removed as duplicate.

3) Reverse-strand soft-clip affects 5′ coordinate (dedup via adjusted 5′):
   - :1143 is chr1 POS=200 CIGAR=45M5S (FLAG 16)
     aligned_end = 200 + 45 - 1 = 244 and 5′ = 244 + 5 = 249
   - :1144 is chr1 POS=200 CIGAR=50M (FLAG 16)
     aligned_end = 200 + 50 - 1 = 249 and  5′ = 249
   Same chr/strand/UMI and same adjusted 5′ so :1144 removed as duplicate.

4) Same adjusted 5′ and UMI but different strand should NOT dedup:
   - :1145 is reverse, chr1 POS=251 CIGAR=50M and 5′ = 300
   - :1146 is forward, chr1 POS=300 CIGAR=50M and 5′ = 300
   Strand is different so both remain.

5) Same chr and POS but different UMI should NOT dedup:
   - :1147 and :1148 at chr1:400 are both kept (UMIs are different).

6) Invalid UMI dropped:
   - :1149 has UMI INVALIDUM so it's removed.

7) Duplicate removal on a second chromosome:
   - :1150 and :1151 are identical on chr2 so only :1150 remains.