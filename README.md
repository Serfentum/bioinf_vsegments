# Finding novel variations of germline immunoglobulin genes using WGS data

This project is dedicated to find new Ig genes from whole genome
sequencing data.


## Methods
Main method of our investigation is aligning. We try to align reads from
WGS data to templates of already known V-genes and find new variants.
Now we are using Biopython built-in local aligner which is similar in
results to aligner from EMBOSS.


## Prerequisites
- python >= 3.6
- Biopython >= 1.70


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
3. Lefranc M-P, Giudicelli V, Duroux P, Jabado-Michaloud J, Folch G, Aouinti S, Carillon E, Duvergey H, Houles A, Paysan-Lafosse T, Hadi-Saljoqi S, Sasorith S, Lefranc G, Kossida S.
IMGT®, the international ImMunoGeneTics information system® 25 years on.
Nucleic Acids Res. 2015 Jan;43(Database issue):D413-22. doi: 10.1093/nar/gku1056. Epub 2014 Nov 5 Free article. PMID: 25378316 LIGM:441


## Acknowledgements
I\`d like to thank everybody especially Bioinformatics Institute.
