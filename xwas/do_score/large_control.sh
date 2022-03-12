

cat ~/athena/ukbiobank/imputed/ukbb.22.sample | cut -f1 -d' ' | tail -n+3 > brit_eid

maxGo=3
rm temp_files/poss_go
for (( i=1; i<=$maxGo; i++ )); do
  echo $i >> temp_files/poss_go
done

echo 1 > temp_files/counter

cat all_specs/author_specs | while read cauthor;do
  low_author=`echo "$cauthor" | tr '[:upper:]' '[:lower:]'`

  for cchr in {1..22};do

    ############################################### SET UP THE GENOTYPE FILE ####################################
    rm temp_files/rsid.${low_author}.${cchr}
    ls ../mod_sets/${low_author} | awk -v var="$cchr" -F'.' '$2 == var {print $0}' | while read cname;do
      cat ../mod_sets/${cauthor}/${cname} | cut -f3 | fgrep rs >> temp_files/rsid.${low_author}.${cchr}
    done

    cat temp_files/rsid.${low_author}.${cchr} | sort | uniq  > temp_files/temp
    mv temp_files/temp temp_files/rsid.${low_author}.${cchr}
    echo start bgenix
    bgenix -g ~/athena/ukbiobank/imputed/ukbb.${cchr}.bgen -incl-rsids temp_files/rsid.${low_author}.${cchr} > temp_files/temp.${low_author}.${cchr}.bgen

    echo start plink2
    plink2_new --memory 12000 --threads 12 --bgen temp_files/temp.${low_author}.${cchr}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${cchr}.sample --keep-fam temp_files/brit_eid --make-bed --out temp_files/geno.${low_author}.${cchr}.almost
    rm temp_files/temp.${low_author}.${cchr}.bgen

    cat temp_files/geno.${low_author}.${cchr}.almost.bim | cut -f2 | sort | uniq -d > temp_files/almost
    plink --bfile temp_files/geno.${low_author}.${cchr}.almost --exclude temp_files/almost --make-bed --out temp_files/geno.${low_author}.${cchr}
    rm temp_files/geno.${low_author}.${cchr}.almost*

    ################################################ DO THE SCORING ###########################################
    ls ../mod_sets/${low_author} | awk -v var="$cchr" -F'.' '$2 == var {print $0}' | while read cname;do
      echo cname $cname
      cchr=`echo $cname | cut -f2 -d'.'`
      cmethod=`echo $cname | cut -f3 -d'.'`
      ver=`echo $cname | cut -f4 -d'.'`

      echo score.${low_author}.${cchr}.${ver}.${cmethod}.profile.zst
      if [ ! -e small_score_files/score.${low_author}.${cchr}.${ver}.${cmethod}.profile.zst ];then

        counter_var=`cat temp_files/counter`
        echo $counter_var
        ./large_score.sh $cname $cauthor $cmethod $cchr $counter_var &> logs/log.${counter_var}.log &
        echo $counter_var + 1 | bc > temp_files/counter


        goOn=False
        while [ $goOn == "False" ]; do
          openSlots=`cat temp_files/poss_go | wc -l`
          sizedir=`du temp_files/ | cut -f1`
          if [ $openSlots -gt 0 ]; then
            if [ $sizedir -lt 50164096 ];then
              echo NOW WE CAN GO
              goOn=True
              sleep $(( ( RANDOM % 3 )  + 1 ))
            fi
          else
            echo MUST WAIT FOR ROOM TO GO
            sleep $(( ( RANDOM % 5 )  + 1 ))
          fi
        done


      else
        echo skipping
      fi

    done

    rm temp_files/geno.${low_author}.${cchr}.*

  done
done
