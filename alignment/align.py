from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from v_segment_generation import combinations



# allele = next(combinations())
# next(combinations('10cores', path_to_heptamers='../data/conserve/hv7',
# path_to_nonamers='../data/conserve/hv9',
# monomers=('A', 'T', 'G', 'C'),
# cum_distribution=(0.25, 0.5, 0.75, 1), length=23))
# reads = []
# with open('/home/arleg/PycharmProjects/igsegments/data/SRR6435693.fasta') as source:
#     readss = SeqIO.parse(source, 'fasta')
#     for i, read in enumerate(readss):
#         reads.append(read)
#         if i > 1000:
#             break
# SeqIO.write(reads, 'some_reads', 'fasta')


allele = 'CTGGGCCTGGACCCAGCAGCCCTCTGGGAAGGCGCTGGGGCACCTCAGCTCCAGGGGCAGCACACACTTCAGCCCAGCCTTTCTGGGCCAACTCTCCATCTGTAGAGACACATCCAAGGCCCAGTTATCCCTGCAGCTGAGCTCCGTGATGGCCAAGGGCAGGGCCGCACATTCCCGTGGGACACAGCG-----------------------ACACAAACG'


def gap_function_for_template(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -100
    return -10


def gap_function_for_read(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -100
    return -10


reads = SeqIO.parse('some_reads', 'fasta')
read1 = next(reads)
rid = read1.description
aa = pairwise2.align.localxc(allele, read1.seq, gap_function_for_template, gap_function_for_read)[0]
for i, j in enumerate(aa):
    print(i, j)




ali1 = MultipleSeqAlignment([SeqRecord(Seq(aa[0], generic_dna), id='v_segment'),
                             SeqRecord(Seq(aa[1], generic_dna), id=rid)])


print(ali1)
print(aa[0])
print(aa[1])
print(rid)
#
# align1 = MultipleSeqAlignment([
#              SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
#              SeqRecord(Seq("ACT-CTAGCTAG", generic_dna), id="Beta"),
#              SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
#          ])


AlignIO.write([ali1], 'al.fa', 'fasta')



