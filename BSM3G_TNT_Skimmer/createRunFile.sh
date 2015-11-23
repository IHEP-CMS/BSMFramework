 #!/bin/bash

folder=$1
datasetName=$2
numberOfJob=$3
proxi=$4

MYWORKAREA=$(pwd)

echo "export X509_USER_PROXY=$MYWORKAREA/../../$proxi"
echo "MYWORKAREA=$CMSSW_BASE/src/"
echo 'cd $MYWORKAREA'
echo 'eval `scram runtime -sh`'
echo "cd $MYWORKAREA"
echo "root -b -q makeSKIM_$numberOfJob.C"

