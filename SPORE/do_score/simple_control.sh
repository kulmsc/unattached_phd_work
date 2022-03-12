
maxGo=4
rm temp_files/poss_go
#cat ../qc/cv_files/train_eid.0.2.txt | cut -f1 > temp_files/brit_eid
#cat ../qc/cv_files/test_eid.0.8.txt | cut -f1 >> temp_files/brit_eid

cat ~/athena/ukbiobank/imputed/ukbb.22.sample | cut -f1 -d' ' | tail -n+3 > brit_eid

for (( i=1; i<=$maxGo; i++ )); do
  echo $i >> temp_files/poss_go
done

echo 1 > temp_files/counter

cat all_specs/author_specs | while read cauthor;do
  cat all_specs/method_specs | while read cmethod;do
    for cchr in {1..22};do

      ls ../mod_sets/${cauthor}/ | fgrep -w ${cmethod} | awk -v var="$cchr" -F. '$2 == var {print $0}' | while read cname;do
        ver=`echo $cname | cut -f4 -d'.'`
        low_author=`echo "$cauthor" | tr '[:upper:]' '[:lower:]'`
        echo score.${low_author}.${cchr}.${ver}.${cmethod}.profile.zst
        if [ ! -e small_score_files/score.${low_author}.${cchr}.${ver}.${cmethod}.profile.zst ];then

          counter_var=`cat temp_files/counter`
          echo $counter_var
          ./simple_score.sh $cname $cauthor $cmethod $cchr $counter_var &> logs/log.${counter_var}.log &
          echo $counter_var + 1 | bc > temp_files/counter

          sleep $(( ( RANDOM % 69 )  + 30 ))

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
              sleep $(( ( RANDOM % 15 )  + 1 ))
            fi
          done
        else
          echo skipping
        fi

      done

    done
  done
done
