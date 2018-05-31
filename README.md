# Finding novel variations of germline immunoglobulin genes using WGS data

This project is dedicated to find new Ig genes from whole genome
sequencing data.


## Methods
Main method of our investigation is aligning. We try to align reads from
WGS data to templates of already known V-genes and find new variants.
Now we are using Biopython built-in local aligner which is similar in
results to aligner from EMBOSS.


## Prerequisites
1. Mandatory - these are necessary to launch
    - python >= 3.6
    - Biopython >= 1.70
2. Facultative - using in simulation
    - numpy >= 1.14.0
    - bowtie2 >= 2.2.6
    - samtools >= 1.7
    - bedtools >= 2.25.0


## Installing
To copy everything into directory **bioinf_vsegments** type in shell:

`git clone https://github.com/Serfentum/bioinf_vsegments.git`


## Usage
For now all fuctionality is inside these scripts:
- align.py - align sequences
- v_segment_generation.py - generate full V segments
- read_generation.py - generate reads with noise from sequence

It is subject to change. We are going to make it available for usage
from cli, so you shouldn\`t try to use until we make everything as it
should be.


## Version
0.1.0

*It is raw version where we are configuring everything on IgH V-segment
to make it work properly.*

*Just wait until release where everything will be done first for IgH
V-segment and afterwards for all other chains.*


## Authors
- Yana Safonova
- Andrey Slabodkin
- Sasha Ilin


## References
1. Greiff, V., Miho, E., Menzel, U., & Reddy, S. T. (2015). Bioinformatic and Statistical Analysis of Adaptive Immune Repertoires. Trends in Immunology, 36(11), 738–749. https://doi.org/10.1016/j.it.2015.09.006
2. Georgiou, G. et. al. (2014). The promise and challenge of high-throughput sequencing of the antibody repertoire, 32(2), 158–168. https://doi.org/10.1038/nbt.2782.The
3. Yaari, G., & Kleinstein, S. H. (2015). Practical guidelines for B-cell receptor repertoire sequencing analysis. Genome Medicine, 7(1), 121. https://doi.org/10.1186/s13073-015-0243-2
1. Lefranc M-P, Giudicelli V, Duroux P, Jabado-Michaloud J, Folch G, Aouinti S, Carillon E, Duvergey H, Houles A, Paysan-Lafosse T, Hadi-Saljoqi S, Sasorith S, Lefranc G, Kossida S. IMGT®, the international ImMunoGeneTics information system® 25 years on. Nucleic Acids Res. 2015 Jan;43(Database issue):D413-22. doi: 10.1093/nar/gku1056. Epub 2014 Nov 5 Free article. PMID: 25378316 LIGM:441
1. Long, H. X., Li, M. Z., & Fu, H. Y. (2016). Determination of optimal parameters of MAFFT program based on BAliBASE3.0 database. SpringerPlus, 5(1). https://doi.org/10.1186/s40064-016-2526-5
1. Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359
1. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
1. Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033
1. Lefranc, M. P., Giudicelli, V., Duroux, P., Jabado-Michaloud, J., Folch, G., Aouinti, S., … Kossida, S. (2015). IMGT R, the international ImMunoGeneTics information system R 25 years on. Nucleic Acids Research, 43(D1), D413–D422. https://doi.org/10.1093/nar/gku1056

## Acknowledgements
I\`d like to thank everybody especially Bioinformatics Institute.
