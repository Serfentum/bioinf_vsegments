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

## alignment
Files with python code for creating whole V-segments and
aligning reads
1. auxiliary - create random fragments, obtain sequences
from given text and fasta files
