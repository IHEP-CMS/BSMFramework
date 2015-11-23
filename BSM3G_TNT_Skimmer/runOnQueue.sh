#!/bin/bash   

folder=$1
queue=$2
proxi=$3

voms-proxy-init -voms cms
cp /tmp/$proxi .
MYWORKAREA=$(pwd)
mkdir $MYWORKAREA/$folder
cd $MYWORKAREA/$folder
cp $MYWORKAREA/fileList.txt .

x=$(cat fileList.txt | wc -l)
for ((i = 1; i < $x+1; i++)) ;
do
    datasetName=($(sed -n "$i"p fileList.txt | awk '{print $1}'))
    rootFilesPath=($(sed -n "$i"p fileList.txt | awk '{print $2}'))
    numerOfFiles=($(sed -n "$i"p fileList.txt | awk '{print $3}'))
    sampleType=($(sed -n "$i"p fileList.txt | awk '{print $4}'))
    echo "Submitting $datasetName sample which has $numerOfFiles files"
    mkdir $MYWORKAREA/$folder/$datasetName
    cd $MYWORKAREA/$folder/$datasetName
    cp $MYWORKAREA/createROOTFile.sh .
    cp $MYWORKAREA/createRunFile.sh .
    chmod 777 createROOTFile.sh
    chmod 777 createRunFile.sh
    for ((j = 1; j < $numerOfFiles+1; j++)) ;
    do
	./createROOTFile.sh $rootFilesPath $datasetName $j $sampleType > makeSKIM_$j.C
	./createRunFile.sh $folder $datasetName $j $proxi > run_$j.sh
	chmod 777 run_$j.sh
	echo "bsub -q $queue run_$j.sh"
	bsub -q $queue run_$j.sh
    done
    cd ..
done