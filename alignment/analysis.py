from Bio import pairwise2
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import re


# Functions are mirrored from package with only one purpose - to make it work on server without issues with
# package path which I have with using __init__.py files. Will fix when will have time
# UP: change this function - add parameter format to make it possible to load fastq files and return SeqRecords
# Place it in alignment directory on ace (of SPAdes))

def get_sequences(path_to_genes='../data/main/simple_fasta/all/ig_hv_all', format='fasta'):
    """
    Read fasta from given path and return generator with all sequences from file
    :param path_to_genes: str - path to fasta file
    :return: generator - generator with sequences
    """
    return SeqIO.parse(path_to_genes, format=format)


def extension(alignment_string, pat=re.compile(r'-+')):
    return max([len(x) for x in re.split(pat, alignment_string)])

path_to_genes = '../data/main/simple_fasta/all/ig_hv_all'
path_to_reads = '../../data/sample.fq'


# Note different formats, cause I've loaded fastq reads
genes = SeqIO.parse(path_to_genes, format='fasta')
reads = get_sequences(path_to_reads, format='fasta')
out = 'alignments.fa'
pat = re.compile(r'-+')
# Yana advise and best result from ROC plot
extension_length = 40
threshold = 20


# Try all combinations of genes and reads to align
# Write appropriate alignments to file
with open(out, 'w') as dest:
    for gene in genes:
        for read in reads:
            for al in pairwise2.align.localxs(read.seq, gene.seq, -5, 0):
                score = al[2]
                s = al[0]
                # Check for alignment score and sequence length
                # Not sure is break for 2nd condition ok and shouldn't be replaced with continue
                if score < threshold or extension(s, pat) < extension_length:
                    break

                # Form alignment and write it to file
                alignment = MultipleSeqAlignment((SeqRecord(Seq(al[1]), id=gene.id, name=gene.name, description=gene.description),
                                                  SeqRecord(Seq(al[0]), id=read.id, name=read.name, description=read.description)),
                                                 alphabet=IUPAC.unambiguous_dna)
                AlignIO.write(alignment, dest, 'fasta')
