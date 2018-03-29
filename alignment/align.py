import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from v_segment_generation import *
from auxiliary.prepare_sequences import *


# According to our score function best align will have score equal to length of template i.e. sequences are identical
# Worst is obviously 0, letters in all positions are different.
def align_one_vs_all(single, other, path_to_out='alignments_of_one_seq_vs_all', threshold=25):
    """
    Align one sequence versus many others, takes SeqRecord objects
    :param single: SeqRecord - template
    :param other: iterable - SeqRecord of reads
    :param path_to_out: str - path to output directory
    :param threshold: float - cutoff whether to align sequences or not
    :return:
    """
    # Preparation - create output directory and compile pattern for editing sequence ids
    os.makedirs(path_to_out, exist_ok=True)
    pattern = re.compile(r'\s+')

    # For each read
    # check score of best alignment, if it greater than cutoff
    # align sequences
    # write alignment to file
    for read in other:
        align_score = pairwise2.align.localms(single.seq, read.seq, 1, -2, -5, -1, score_only=True)

        if align_score >= threshold:
            align = pairwise2.align.localms(single.seq, read.seq, 1, -2, -5, -1)[0]
            als = MultipleSeqAlignment([SeqRecord(Seq(align[0], generic_dna), id=single.id),
                                        SeqRecord(Seq(align[1], generic_dna), id=read.id)])

            AlignIO.write([als], f'{path_to_out}/{re.sub(pattern, "_", read.id)}-{align_score}', 'fasta')


def align_all_vs_all(seqs_a, seqs_b, path_to_out='alignments_of_all_vs_all', threshold=25):
    """
    Align each sequence from seqs_a versus all sequences from seqs_b, takes SeqRecord objects.
    Write to files pairwise alignments with score greater or equal than cutoff.
    :param seqs_a: iterable - pack of SeqRecord objects, by convention let`s supply here templates
    :param seqs_b: iterable - pack of SeqRecord objects, by convention let`s supply here reads
    :param path_to_out: str - path to output directory
    :param threshold: float - cutoff whether to align sequences or not
    :return:
    """
    # Create directory and compile pattern for naming subdirectories
    os.makedirs(path_to_out, exist_ok=True)
    pattern = re.compile(r'\s+')

    for seq in seqs_a:
        align_one_vs_all(seq, seqs_b, f'{path_to_out}/{re.sub(pattern, "_", seq.id)}', threshold)


# 1st try with just 1 template made by hand
v_segment = 'CTGGGCCTGGACCCAGCAGCCCTCTGGGAAGGCGCTGGGGCACCTCAGCTCCAGGGGCAGCACACACTTCAGCCCAGCCTTTCTGGGCCAACTCTCCATCTGTAGAGACACATCCAAGGCCCAGTTATCCCTGCAGCTGAGCTCCGTGATGGCCAAGGGCAGGGCCGCACATTCCCGTGGGACACAGCG-----------------------ACACAAACG'
single = SeqRecord(Seq(v_segment), id='template')
# Read few reads
reads = SeqIO.parse('some_reads', 'fasta')
# Read 10 v_segments
vss, reads = prepare_sequences('10cores', 'some_reads')

align_all_vs_all(vss, reads, 'ala', 12)


# allele = 'CTGGGCCTGGACCCAGCAGCCCTCTGGGAAGGCGCTGGGGCACCTCAGCTCCAGGGGCAGCACACACTTCAGCCCAGCCTTTCTGGGCCAACTCTCCATCTGTAGAGACACATCCAAGGCCCAGTTATCCCTGCAGCTGAGCTCCGTGATGGCCAAGGGCAGGGCCGCACATTCCCGTGGGACACAGCG-----------------------ACACAAACG'
# s = 'AATAACATTGATACTACATACCATGGTTTCACTGCATATGAAAAAATAAAAGATGATTTGTTCTAACTTTAAACATATGCACTTTCTGTTGATCTACTGTACCTCAATAGAACTGTTTTAAAATAAAAATTACAAAATTATAAGATTTATAGGTTTTAAGGTTTTATCACAGAGCAGATTTACCATAAGAAACCACAATTTCCCAAATGCTATCAATATCACAAATCTCCCCAGGACACTGTCACGTGCTCTGAGCCCCACTCTCTCCAAAGGCCTCTAACCAGAGAGCTTACTATATAGTAGGAGACATGGAAATAGAGCCCTCCCTCTGCTTATGAAAACCAGCCCAGCCCTGACCCTGCAGCTCTGGGACAGGAGCCCCAGCCCTGGGATTTTCAGGTGTTTTCATTTGGTGATCAGGACTGAACACAGAGGACCACCAAGGAGTCATGGCTGAGCTGGCTTTTTCTTGTGGCTATTTTAAAAGGTAATTCATGGAGAAATAGAAAAATTGAGTGTGAGTGGATAAGAGTGAGATAAACAGTGGATTTGTGTGGAAGTTTCTGACCAGGTTGTCTCTTTGTTTGCAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTAAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTGACTACTACATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTACCATATACTACGCAGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGACACAGTGAGGGGAAGTCAGTGTGAGCCCAGACACAAACCTCCCTGCAGGGGTCCCCAGGACCACCAGGGGGCGCCCGGGACACTGTGCACGGGGCTGTCTCCAGGGCAGGTGCAGGTGCTGCTGAGGCCTGGCTTCCCTGTCATGGCCTGGGCGGCCTCGTTGTCAAATTTCTCCAGGGAACTTCTCCAGATTTACAATTCTGTACTGACATTTCATGTCTCTAAATGCAAAACTTTTTTGTTCTTTTTGTATTTTTGTTTTTGTAACAGGAGGACACACCCTCACCTCCACAGAAGCCACAGTGTCACTTTGGGGGCAGAT'




