chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

ess=`cat temp_files/ss.${low_author}.${chr} | head -2 | tail -n +2 | cut -f9`
eth=`cat helper_scripts/author_conversion | fgrep $author | cut -f2 -d' '`

plink --bfile ~/athena/refs/1000genomes/${eth}.${chr} --freq --out temp_files/${low_author}.${chr}

Rscript helper_scripts/add_maf.R $d $eth $author $chr


python ~/Programs/polyfun/munge_polyfun_sumstats.py \
    --sumstats ${d}/ss.${low_author}.${chr}.maf \
    --out ${d}/munged.ss \
    --min-info 0 \
    --n $ess


#create fine-mapping jobs
python ~/Programs/polyfun/create_finemapper_jobs.py \
    --sumstats ${d}/munged.ss \
    --n $ess \
    --geno ~/athena/refs/1000genomes/${eth}.${chr} \
    --method susie \
    --non-funct \
    --allow-missing \
    --max-num-causal 5 \
    --out-prefix ${d}/polyfun_output \
    --jobs-file ${d}/jobs.txt 


bash ${d}/jobs.txt


#aggregate all of the results
python ~/Programs/polyfun/aggregate_finemapper_results.py \
    --out-prefix ${d}/polyfun_output \
    --sumstats ${d}/munged.ss \
    --out ${d}/polyfun_output.agg.txt.gz \
    --adjust-beta-freq



Rscript helper_scripts/polypred_add_beta.R $d $eth $author $chr 1

for curr_eth in european hispanic african eastasian;do
  Rscript helper_scripts/combo_polypred_prscs.R $chr $curr_eth 0 2
  Rscript helper_scripts/combo_polypred_prscs.R $chr $curr_eth 3 4
  Rscript helper_scripts/combo_polypred_prscs.R $chr $curr_eth 6 6
done

Rscript helper_scripts/combo_polypred_prscs.R $chr total 1 2
Rscript helper_scripts/combo_polypred_prscs.R $chr total 4 4
Rscript helper_scripts/combo_polypred_prscs.R $chr total 7 6

 

#later on need to combine linearly with prscs computations
