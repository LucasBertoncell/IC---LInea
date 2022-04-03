#!/bin/bash

# To give a bash file permission of execution: chmod u+x filename.sh

##  BASH ##


echo Hello. This is your first bash script

# JACKKNIFE SUB-SAMPLES ARE HERE: /jackk_cats
# SAVE TO: /CUTE/sub_samples

end=$((SECONDS+3))

while [ $SECONDS -lt $end ]; do
    :/nohup condor_do CUTE dinamic_ini.c 

done