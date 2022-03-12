author=Bentham
sec_author=Christophersen
low_author=bentham
o_author=christophersen
Rscript ~/Programs/HDL/HDL.run.R \
    gwas1.df=${author}/${low_author}.munged.sumstats.gz \
    gwas.2.df=${sec_author}/${o_author}.munged.sumstats.gz
    LD.path=~/athena/refs/UKB_imputed_SVD_eigen99_extraction \
    output.file=gen_corr/${low_author}.${o_author}.corr.hdl

