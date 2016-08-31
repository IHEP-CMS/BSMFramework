#!/bin/bash
#list=(First2016_TTHbb First2016_TTHnbb First2016_TT First2016_2_ST First2016_2_SaT First2016_2_DY First2016_2_WJets First2016_2_WW First2016_2_WZ First2016_2_ZZ)
list=(Firts2016_SEleB Firts2016_SEleC First2016_SMuB First2016_SMuC)
pos=0
for d in ${list[@]}; 
do
  #echo "${list[$pos]}"
  crab status $d/crab_$d
  let pos=pos+1
done
