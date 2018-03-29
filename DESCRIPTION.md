# Description
Information about stored files.

## data
Files with data about genes etc.

Pattern for files - ig_*v or *v where star is a
placeholder for name of chain.

### Structure
1. conserve - folder for conserved sequences - heptamers
and nonamers

    - hv7 - IgH V-segment heptamer
    - hv9 - IgH V-segment nonamer
    - And so on

2. main - folder for cores of V segment gene
    1. aligned - folder for loaded from IMGT alignment of
    core genes

        - ig_hv - alignment of IgH V-segment genes
        - ig_hklv - alignment of IgH, K and L V-segment
        genes
        - And so on
    2. simple_fasta - folder for fasta sequences of
    Ig V segment genes

        1. all - all genes including out of ORF
        pseudogenes
            - ig_hv_all - all IgH V-segment genes
            - And so on
        2. without_outframe_pseudogenes - all genes
        except out of ORF pseudogenes
            - ig_hv_wop - IgH V-segment genes without
            out of ORF pseudogenes
            - And so on
    3. index - folder with files of IgH locus index created by bowtie2
    4. intervals - folder with files from alignment V genes on reference
    5. reference - folder with IgH locus reference

## alignment
Files with python code for creating whole V-segments and
aligning reads
1. v_segment_generation.py - create full V segments
2. align.py - align sequences, obtain SeqRecords from fasta

## auxiliary
Files with python code for everything supportive
1. create_spacer.py - create random fragments
2. get_sequences.py - obtain sequences from fasta
3. read_sequences.py - obtain sequences from text file
4. heptamers.py - obtain heptamers from text file
5. nonamers.py - obtain nonamers from text file
6. prepare_sequences.py - read 2 fasta
7. sequencing.py - sequence 1 read from sequence or it whole

## preprocessing
Scripts that should be launched as you get data
1. sequence_capitilization.py - shift all sequence in fasta file to
upper case
2. indexing.py - index V genes to Ig locus and get file with gene
intervals

## simulation
read_generation.py - create reads from sequence separate reads
overlapping with V genes