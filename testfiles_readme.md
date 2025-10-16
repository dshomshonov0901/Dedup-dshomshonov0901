The first two reads are identical across chromosome, strand, position, and UMI.

The third read has same UMI but reverse strand (FLAG 16) and different CIGAR (soft clipping), producing a different adjusted 5′ start → retained.

The fourth and fifth reads have different positions/UMIs → retained.

The last read has an invalid UMI → discarded.
