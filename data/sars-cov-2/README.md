# SARS-CoV-2 data

- `reference.fa`: Wuhan-Hu-1 ([MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3/))
- `reference.gff`: Modified version of [GCF_009858895.2_ASM985889v3_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz). Replaced `NC_045512.2` with `MN908947.3` (identical entries in NCBI; RefSeq vs GenBank). Using the RefSeq GFF (NC_045512.2) instead of the version that one can download from NCBI Nucleotide DB ( GenBank MN908947.3), since the RefSeq GFF has more entries and info, which may be useful for variant effect quantification.
- `primer_schemes/Freed_2020.bed`: Downloaded from archive listed in [Freed et al (2020)](https://academic.oup.com/biomethods/article/5/1/bpaa014/5873518) [Zenodo](https://zenodo.org/record/3897530)
- `primer_schemes/ARTIC_V3.bed`: Dowloaded from [artic-network/artic-ncov2019](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) GitHub (commit [87449a03f5d85230a2a9214064adf2e2cc2abd51](https://github.com/artic-network/artic-ncov2019/commit/87449a03f5d85230a2a9214064adf2e2cc2abd51))
