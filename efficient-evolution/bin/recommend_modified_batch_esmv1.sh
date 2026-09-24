#/bin/bash
# usage: efficient-evolution batch Script
# ./recommend_modified_esmv1.sh ~/inputs/
# ~/inputs/ folder contain fasta files, suffix of files can be .fa or .fasta
# alpha = 1 and K = 1 by default
input_folder=$1
echo "processing fasta files in $input_folder"
# modify python and recommend_modified.py paths if needed
ls $input_folder | grep -P "fa$|fasta" | while read input_file
do
	bioawk -c fastx '{cmd="/home/zhangxiang/anaconda3/envs/msa-transformer/bin/python3 /home/zhangxiang/gitFolder/CasDelta/efficient-evolution/bin/recommend_modified.py --gene "$name" --sequence "$seq" --model-version esmv1 --topK 1";system(cmd);}' ${input_folder}/${input_file}
done



