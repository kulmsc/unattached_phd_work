chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


if [ $author == european ];then
  eth=eur
elif [ $author == african ];then
  eth=afr
elif [ $author == eastasian ];then
  eth=eas
elif [ $author == hispanic ];then
  eth=hsp
elif [ $author == total ];then
  eth=tot
fi

Rscript helper_scripts/get_impact_snps.R $eth $d $chr

i=1
cat all_specs/impact_param_specs | tail -n +2 | while read spec;do
  plim=`echo $spec | cut -f1 -d' '`
  r2lim=`echo $spec | cut -f2 -d' '`
  perclim=`echo $spec | cut -f3 -d' '`

  #if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.clump.${i}.ss ]; then
    plink --memory 4000 --threads 1 --bfile geno_files/${low_author}.${chr} --extract ${d}/rsid.${perclim}.txt --make-bed --out temp_files/${low_author}.${chr}.small

    plink --memory 4000 --threads 1 --bfile temp_files/${low_author}.${chr}.small --clump temp_files/ss.${low_author}.${chr} --clump-snp-field RSID --clump-p1 $plim --clump-r2 $r2lim --out ${d}/out

    if [ -f ${d}/out.clumped ]; then
      sed -e 's/ [ ]*/\t/g' ${d}/out.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > ${d}/done_rsids
      fgrep -w -f ${d}/done_rsids temp_files/ss.${low_author}.${chr} > ~/athena/SPORE/mod_sets/${author}/${low_author}.${chr}.impact.${i}.ss
    fi
  #fi

  rm temp_files/${low_author}.${chr}.small.bed temp_files/${low_author}.${chr}.small.bim temp_files/${low_author}.${chr}.small.fam

  let i=i+1
done





#go through all annotation cols and determine the one with the largest Prop_h2 whose Prop_h2 - Prop_h2_se is greater than zero
#pull that column from that matrix of impact scores
#-- need to repeat for each ancestry --
#then do clumping for a predefined set of SNPs who are in a top percentile of the impact score
#typically do top 5% IMPACT SNPs, can do a range along with normal clumping parameters

#-------------------------------------------

#FROM THE PAPER
# For functionally informed PRS, we
#restricted the analysis to variants with high IMPACT score according to the lead IMPACT
#annotation before conducting LD clumping. As before, we define the lead annotation as the
#one with the largest τ* estimate that was significantly greater than 0. When we designed
#PRS-EUR, we utilized the lead IMPACT annotation in EUR GWAS summary statistics
#
#We applied S-LDSC (v.1.0.0)3
# to partition the
#common (MAF > 5%) SNP heritability of 111 polygenic traits and diseases. We partitioned
#heritability using a customized version of the baseline-LD model, accounting for 69 celltype-nonspecific baseline-LD annotations, and added one or more IMPACT annotations to
#the model to test for cell-type-specific enrichment. 
