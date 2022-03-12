
maxDirs=8


	rm temp_files/*
	rm logs/*
	rm done_check/*
	rm -r comp_zone/*

	for (( i=1; i<=$maxDirs; i++ )); do
		mkdir comp_zone/dir$i
		echo $i >> temp_files/poss_dirs
	done



	cat all_specs/chr_specs | while read chr; do

		cat all_specs/author_specs | while read author; do		

			./setup_data.sh $author $chr

			cat all_specs/method_specs | while read method; do


				echo about to hermes
				./hermes.sh $chr $method $author &> logs/${author}.${chr}.${method}.log &
				sleep $(( ( RANDOM % 30 )  + 1 ))

				goOn=False
				while [ $goOn == "False" ]; do
					openSlots=`cat temp_files/poss_dirs | wc -l`
					if [ $openSlots -gt 0 ]; then
						echo NOW WE CAN GO
						goOn=True
					else
						echo MUST WAIT FOR ROOM TO GO
						sleep $(( ( RANDOM % 30 )  + 1 ))
						openSlots=`cat temp_files/poss_dirs | wc -l`
					fi
				done

			done		

			./check_to_delete.sh

		done
	done
