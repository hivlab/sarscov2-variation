#!/bin/bash
fileid="0B3llHR93L14wd0pSSnFULUlhcUk"
filename="hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
