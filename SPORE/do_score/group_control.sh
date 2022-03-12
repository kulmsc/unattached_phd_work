#! /bin/bash -l

#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=rgwas
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T

spack load -r /mjrrusu

echo 1 > temp_files/counter

cat all_specs/author_specs | while read cauthor;do
  cat all_specs/method_specs | while read cmethod;do
    for cchr in {5..22};do

      rm temp_files/all_rsid temp_files/all_modsets
      ls ../mod_sets/${cauthor}/ | fgrep -w ${cmethod} | awk -v var="$cchr" -F. '$2 == var {print $0}' | while read cname;do
        if [ ! -e small_score_files/score.${low_author}.${cchr}.${ver}.${cmethod}.profile.zst ];then
          echo $cname >> temp_files/all_modsets
          cat ../mod_sets/${cauthor}/${cname} | tail -n+2 | cut -f3 >> temp_files/all_rsid
        fi
      done

      cat temp_files/all_rsid | sort | uniq > temp_files/temp; mv temp_files/temp temp_files/all_rsid
      bgenix -g ~/athena/ukbiobank/imputed/ukbb.${cchr}.bgen -incl-rsids temp_files/all_rsid > temp_files/temp.bgen

      Rscript combine_modsets.R

      num_cols=`head -1 temp_files/comboset.txt | cut -f3-300 | tr '\t' '\n' | wc -l`

      plink2 --bgen temp_files/temp.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${cchr}.sample --score temp_files/comboset.txt 1 2 header-read cols=+scoresums --score-col-nums 4-${num_cols}  --out combo_scores/${cauthor}.${cmethod}.${cchr}

      rm temp_files/comboset.txt temp_files/all_rsid temp_files/temp.bgen temp_files/all_modsets

    done
  done
done
