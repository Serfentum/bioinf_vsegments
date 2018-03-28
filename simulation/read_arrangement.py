import py


# Coordinates for IgH locus was taken according to this article and quite expanded
# https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=3492
# Left end a - 105392425
# Right end b - 107043718

# After bowtie2-samtools-bedtools coordinate numeration starts from 0 and going from right to left
# Because of this for example gene IGHV(III)-82*01
# https://www.ncbi.nlm.nih.gov/gene/28338
# which has coordinates 106879553-106879844 in reference have been shifted to 163874-164166
# Thats why we make this rearrangement



