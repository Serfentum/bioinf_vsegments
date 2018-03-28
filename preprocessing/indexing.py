import os
import subprocess as sp
import shlex


def gene_alignment_bt2_sam_bed(genes='../data/main/simple_fasta/all/ig_hv_all',
                               reference='../data/reference/igh_locus_minus_strand',
                               output_path='preintervals'):
    """
    Obtain reference intervals where genes dwells. Create folder with index, sam, bam, bed and preintervals file with
    genes starts and ends.
    :param genes: str - path to aligning genes
    :param reference: str - path to reference genome
    :param output_path: str - path to dump everything
    :return:
    """
    # Count available cores and create output directories
    cores_number = len(os.sched_getaffinity(0))
    os.makedirs(f'{output_path}/reference_index', exist_ok=True)

    # Commands to do
    # indexing reference genome -> directory with bt2 files
    # aligning core genes       -> sam file
    # converting sam to bam     -> bam file
    # preparing sorted bam      -> sorted bam file
    commands = (f'bowtie2-build {reference} {output_path}/reference_index/ref -p {cores_number}',
                f'bowtie2 -f -x {output_path}/reference_index/ref -U {genes} -S {output_path}/genes.sam -p {cores_number}',
                f'samtools view -b {output_path}/genes.sam -o {output_path}/genes.bam -@ {cores_number}',
                f'samtools sort {output_path}/genes.bam -o {output_path}/genes.sorted.bam -@ {cores_number}')

    # Prepare commands
    commands = map(lambda command: shlex.split(command), commands)
    # Execute commands
    for command in commands:
        sp.check_call(command, universal_newlines=True)

    # Create bed file from bam
    # Write name of reference, start and end intervals to tsv file
    command_bedtools = shlex.split(f'bedtools bamtobed -i {output_path}/genes.sorted.bam')
    with open(f'{output_path}/genes.bed', 'a+') as bed, open(f'{output_path}/preintervals', 'w') as dest:
        sp.check_call(command_bedtools, universal_newlines=True, stdout=bed)
        bed.seek(0)
        for line in bed:
            line = line.split('\t')
            dest.write(f'{line[3]}\t{line[1]}\t{line[2]}\n')

