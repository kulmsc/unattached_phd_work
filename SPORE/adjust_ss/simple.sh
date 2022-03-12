chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


./setup_data.sh $author $chr
./ldak.sh $chr $author $dir
rm geno_files/*
rm ${d}/*




#if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.simple.1.ss ]; then
#  cat temp_files/ss.${low_author}.${chr} | awk '$8 < 5e-8 {print $0}' > ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.simple.1.ss
#fi

#if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.simple.2.ss ]; then
#  cat temp_files/ss.${low_author}.${chr} | awk '$8 < 1e-6 {print $0}' > ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.simple.2.ss
#fi

#if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.simple.3.ss ]; then
#  cat temp_files/ss.${low_author}.${chr} > ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.simple.3.ss
#fi

#if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.simple.4.ss ]; then
#  cat temp_files/ss.${low_author}.${chr} | awk '$8 > 0.5 {print $0}' > ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.simple.4.ss
#fi
