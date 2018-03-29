import re
from auxiliary.sequencing import sequence_read
from auxiliary.get_sequences import get_sequences
from auxiliary.read_sequences import read_sequences


# I don`t know why but I have a plethora of leading Ns in a reference sequence
# Obviously reads from them won`t be very diverse, so it`d be better to skip them
# Yet we need indices of nucleotides be consistent with numbering in our file with intervals


def sequencing_simulation(intervals='/home/arleg/ig_construction/data/main/intervals/intervals',
                          sequence='/home/arleg/ig_construction/data/reference/igh_locus_minus_strand.fa',
                          interlap_length=10, read_length=100, error_rate=0.1):
    """
    Simulate sequencing upon given sequence with specified read length and error rate.
    Reads which overlaps with genes from intervals are marked and return separately from non-overlapping.
    Writes file with non-overlapped unique reads and unique overlapped with genes reads.
    :param intervals: str - path to tsv file with 3 columns - name, start and end of gene in coordinates of given sequence
    :param sequence: str - path to fasta with sequence, Ig in our case
    :param interlap_length: int - number of intersected bases in read and gene to treat read belonging to gene
    :param read_length: int - length of reads
    :param error_rate: float - fraction of errors
    :return:
    """
    # Read intervals data and extract very intervals
    intervals = (interval.split('\t') for interval in read_sequences(intervals))
    intervals = sorted(list({(int(start), int(stop)) for _, start, stop in intervals}), key=lambda x: x[0])
    # Read sequence
    sequence = next(get_sequences(sequence))

    # Find 1st sense nucleotide in it
    pattern = re.compile(r'[ATGC]')
    start = re.search(pattern, sequence).start()

    # Create sets for reads
    overlapped = set()
    non_overlapped = set()

    # Iterate over sequence
    # look at intervals
    # if read will be overlapped with gene, sequence all overlapped with this gene reads and add them to both sets
    # otherwise sequence read and add to all reads
    it = iter(range(start, len(sequence) - read_length + 1))
    for i in it:
        for interval in intervals:
            while interval[0] + interlap_length <= i < interval[1] - interlap_length:
                read = sequence_read(sequence, i, read_length, error_rate)
                overlapped.add(read)
                i = next(it)
        read = sequence_read(sequence, i, read_length, error_rate)
        non_overlapped.add(read)

    return overlapped, non_overlapped


# a = sequencing_simulation('/home/arleg/ig_construction/alignment/test.int',
#                           '/home/arleg/ig_construction/alignment/test.fa',
#                           2, 4, 0)

a = sequencing_simulation()

with open('overlapped', 'w') as dest:
    for read in a[0]:
        dest.write(f'{read}\n')
print(f'\n***\n')
with open('non_overlapped', 'w') as dest:
    for read in a[1]:
        dest.write(f'{read}\n')


