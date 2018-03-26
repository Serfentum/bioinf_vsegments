from Bio import SeqIO


def sequence_capitalization(path):
    """
    Rewrite sequences in fasta to make them all in upper case
    :param path: str - path to file
    :return:
    """
    seqs = list(SeqIO.parse(path, 'fasta'))
    for seq in seqs:
        seq.seq = seq.seq.upper()
    SeqIO.write(seqs, path, 'fasta')

# sequence_capitalization('/home/arleg/ig_construction/data/main/simple_fasta/all/ig_hv_all')
