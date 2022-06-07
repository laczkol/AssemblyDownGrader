#!/bin/bash

set -euo pipefail -euo nounset

infile=" "
minpoints="1"
maxpoints="5"
minrep="1"
maxrep="3"
mingap="1"
maxgap="20"
minlen="1"
outprefix="draft"
sn="0.0002"
in="0.0001"
de="0.0001"
iv="0.00005"
inmax="10"
ivmax="100"
demax="5"
du="0"
tl="0"
dumin="1"
dumax="2"
tlmin="1"
tlmax="2"

while [[ "$#" -gt 0 ]];
do
    case $1 in
        -i|--infile)
            infile="$2"
            shift
            ;;
        -minp|--minpoints)
            minpoints="$2"
            shift
            ;;
        -maxp|--maxpoints)
            maxpoints="$2"
            shift
            ;;
        -r|--reps)
            maxrep="$2"
            shift
            ;;
        -ming|--mingaps)
            mingap="$2"
            shift
            ;;
        -maxg|--maxgaps)
            maxgap="$2"
            shift
            ;;
        -minlen|--minimum-length)
            minlen="$2"
            shift
            ;;
        -out|--out-prefix)
            outprefix="$2"
            shift
            ;;
        -sn|--snp-rate)
            sn="$2"
            shift
            ;;
        -in|--insertion-rate)
            in="$2"
            shift
            ;;
        -de|--deletion-rate)
            de="$2"
            shift
            ;;
        -iv|--inversion-rate)
            iv="$2"
            shift
            ;;
        --inmax|--insert-max-length)
            inmax="$2"
            shift
            ;;
        --ivmax|--inversion-max-length)
            ivmax="$2"
            shift
            ;;
        --demax|--deletion-max-length)
            demax="$2"
            shift
            ;;
        --du|--duplication-rate)
            du="$2"
            shift
            ;;
        --tl|--translocation-rate)
            tl="$2"
            shift
            ;;
        --dumin|--duplication-min-length)
            dumin="$2"
            shift
            ;;
        --dumax|--duplication-max-length)
            dumax="$2"
            shift
            ;;
        --tlmin|--translocation-min-length)
            tlmin="$2"
            shift
            ;;
        --tlmax|--translocation-max-length)
            tlmax="$2"
            shift
            ;;
        *) echo "Unknown parameter passed: $1"
           exit 1
           ;;
    esac
    shift
done

if [[ ${#infile} -le 1 ]]; then
    echo "Input file is not specified or empty"
    exit 1
fi

cat <<'END_FIGLET'
    _                           _     _       ____                       ____               _
   / \   ___ ___  ___ _ __ ___ | |__ | |_   _|  _ \  _____      ___ __  / ___|_ __ __ _  __| | ___ _ __
  / _ \ / __/ __|/ _ \ '_ ` _ \| '_ \| | | | | | | |/ _ \ \ /\ / / '_ \| |  _| '__/ _` |/ _` |/ _ \ '__|
 / ___ \\__ \__ \  __/ | | | | | |_) | | |_| | |_| | (_) \ V  V /| | | | |_| | | | (_| | (_| |  __/ |
/_/   \_\___/___/\___|_| |_| |_|_.__/|_|\__, |____/ \___/ \_/\_/ |_| |_|\____|_|  \__,_|\__,_|\___|_|
                                        |___/                                                       v0.2
END_FIGLET

echo "Run started at $(date)"

rep=$(seq "$minrep" "$maxrep")

bioawk -c fastx '{ print $name }' "$infile" > names

sname=$(cut -f 1 names)
fname=$(echo "$infile" | cut -f 1 -d ".")

for r in $rep
do

    echo " "
    echo "#################"
    echo "#### Repat $r ####"
    echo "#################"
    echo " "

    cp "${infile}" r"${r}"_"${fname}".fa

    mutation-simulator r"${r}"_"${fname}".fa args -sn "$sn"\
        -in "$in"\
        -de "$de"\
        -iv "$iv"\
        -inmax "$inmax"\
        -ivmax "$ivmax"\
        -demax "$demax"\
        -du "$du"\
        -tl "$tl"\
        -dumin "$dumin"\
        -dumax "$dumax"\
        -tlmin "$tlmin"\
        -tlmax "$tlmax"

    sed -i 's/>>/>/g' r"${r}"_"${fname}"_ms.fa #to correct double headers sometimes added by mutation simulator, strange

    bioawk -c fastx '{ print $name, length($seq) }' r"${r}"_"${fname}"_ms.fa > name_length

    for i in $sname
    do

        echo "#### Breaking up contig $i ####"

        slen=$(grep -w "$i" name_length | cut -f 2)
        numpoints=$(shuf -i "${minpoints}"-"${maxpoints}" -n 1)

        start="0" #start coordinate, can be hardcoded

        shuf -i "${start}"-"${slen}" -n "$numpoints" | sed -e "1i${start}" -e "1i${slen}" | sort -n  > numbers

        gapnum=$(shuf -i "$mingap"-"$maxgap" -n 1)

        paste <(awk -v slen="$slen" '$1 < slen' numbers) <(awk '$1 > 0' numbers) |\
            perl -pe "s/^/${i}\t/" |\
            while read line
            do
                echo "$line" |\
                    awk -v OFS="\t" -v gap="$gapnum" '$3 = $3 - gap'
            done |\
            awk -v OFS="\t" '$4 = $3 - $2' |\
            while read line
            do
                echo "$line" |\
                    #creates overlaps; alternative is to re-add $gapnum awk -v OFS="\t" -v gap="$gapnum" '{if ($4 < 0) $3 = $3 + gap}'
                    awk -v OFS="\t" -v slen="$slen" '{if ($4 < 0) $3 = $3 - $4*2; $4 = $3 - $2; if ($3 > slen) $3 = slen; print $0}' |\
                        awk -v OFS="\t" -v minlen="$minlen" '$4 > minlen'
            done  > temp_"${i}"_r"${r}".bed

    done

    #uniqid=$(date | sed -e 's/ /_/g' -e 's/\.//g' -e 's/,//g' -e 's/://g' | cut -f 1,2,3,5 -d_)

    uniqid=$(date | tr -cd [:alnum:])

    cat temp*_r"${r}".bed > "${infile}"_r"${r}"_"${uniqid}".bed

    echo " "
    echo "#### Writing contig sequences of $infile to ${outprefix}_r${r}_${uniqid}_${infile} ####"

    bedtools getfasta -fi r"${r}"_"${fname}"_ms.fa\
        -bed "${infile}"_r"${r}"_"${uniqid}".bed > "${outprefix}"_r"${r}"_"${uniqid}"_"${infile}"

done

rm temp* #to tidy up
rm *fa.fai #so the next run does not throw errors

echo "Run ended at $(date)"


#TO DO
#scale number of breaks with length
