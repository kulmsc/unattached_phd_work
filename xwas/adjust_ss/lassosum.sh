chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}



echo started shell script

Rscript helper_scripts/lassosum.R $d $chr $low_author


