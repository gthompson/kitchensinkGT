#!/bin/bash
cd $1
for yyyy in ????/
do
    cd $yyyy
    for mm in ??/
        do
            echo $yyyy $mm `ls ${mm}/*.mseed | wc -l` `ls -l ${mm}/*.csv | wc -l`
        done
    cd ..
done
