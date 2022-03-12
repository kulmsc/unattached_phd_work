#! /bin/bash -l

#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=rgwas
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=16G   # memory requested, units available: K,M,G,T




eth=hsp
fulleth=hispanic

for impact_col in {10..12};do
#impact_col=$1
let real_col=impact_col-4


for chr in {1..22};do

  rm temp temp2 temp3 curr_annot current.bim current.bed current.fam curr_rsids

  echo 1 CHR: $chr

  zcat ready_files/${eth}.${chr}.annot.gz | tail -n+2 | cut -f3 > other_annot_rsids
  zcat ~/athena/exome_score/funct_anno/impact/IMPACT707_EUR_chr${chr}.annot.gz | cut -f1-4,${impact_col} | fgrep -w -f other_annot_rsids > curr_annot
  cat curr_annot | cut -f3 > curr_rsids

  echo 2 CHR: $chr

  plink --bfile ~/athena/refs/1000genomes/${eth}.${chr} --extract curr_rsids --make-bed --out current
  cat current.bim | wc -l
  cat current.bim | cut -f2 > temp
  zcat ~/athena/exome_score/funct_anno/impact/IMPACT707_EUR_chr${chr}.annot.gz | cut -f1-4,${impact_col} | fgrep -w -f temp > curr_annot

  echo 3 CHR: $chr


  cat curr_annot | cut -f1 > temp
  cat curr_annot | cut -f3-100 > temp2
  cat current.bim | cut -f4 > temp3
  echo "CHR BP SNP CM base" | tr ' ' '\t' > curr_annot #the col heads should be CHR	BP	SNP	CM	base
  paste temp temp3 temp2 >> curr_annot

  echo 4 CHR: $chr
  cat current.bim | wc -l
  head curr_annot

  taskset -c 1-4 ldsc\
   --l2\
   --bfile current\
   --ld-wind-cm 1\
   --annot curr_annot\
   --out impact_files/curr_impact.${eth}.${chr}

  #do not need to manually move impact_files/curr_impact* because of ln -s
  mv curr_annot ready_files/curr_impact.${eth}.${chr}.annot
  gzip -f ready_files/curr_impact.${eth}.${chr}.annot

  mv impact_files/curr_impact.${eth}.${chr}.l2.ldscore.gz ready_files/curr_impact.${eth}.${chr}.l2.ldscore.gz

  rm temp temp2 temp3 curr_annot current.bim current.bed current.fam curr_rsids

  echo CHR: $chr
  echo BOTTOM

done

cd ready_files
for i in {1..22};do
  rm curr_impact.${eth}.${i}.l2.M
  rm curr_impact.${eth}.${i}.l2.M_5_50
  rm curr_impact.${eth}.${i}.l2.ldscore.gz

  ln -s ../impact_files/curr_impact.${eth}.${i}.l2.M .
  ln -s ../impact_files/curr_impact.${eth}.${i}.l2.M_5_50 .
  ln -s ../impact_files/curr_impact.${eth}.${i}.l2.ldscore.gz .

done
cd ..

taskset -c 1-4 ldsc \
  --h2 ~/athena/SPORE/raw_ss/${fulleth}/${fulleth}.munged.sumstats.gz \
  --ref-ld-chr ready_files/${eth}.,ready_files/curr_impact.${eth}. \
  --frqfile-chr freq_files/${eth}. \
  --w-ld-chr weight_files/w.${eth}. \
  --overlap-annot \
  --print-coefficients \
  --print-delete-vals \
  --out BMI.baselineLD_yourannot

mv BMI.baselineLD_yourannot.results final_results/annot.${eth}.${real_col}
rm BMI.baseline*

###

#taskset -c 40-44 ldsc --h2 ~/athena/SPORE/raw_ss/european/european.munged.sumstats.gz --ref-ld-chr ready_files/${eth}.,ready_files/curr_impact.${eth}. --frqfile-chr freq_files/${eth}. --w-ld-chr weight_files/w.${eth}. --overlap-annot --print-coefficients --print-delete-vals --out BMI.baselineLD_yourannot

done
