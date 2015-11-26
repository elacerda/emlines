#!/bin/bash
#dryrun=true
dryrun=false
#for file in /Users/lacerda/CALIFA/listOf298GalPrefixes.txt
#for file in debug.txt
for file in /Users/lacerda/CALIFA/listv20_q050.d15a.txt
#for file in /Users/lacerda/CALIFA/list_S[abcd]*.txt
do
    bname=$(basename $file)
    name=${bname%.txt}
    global_opts="--nolinecuts"
    rbinstep=0.1
    
    ##########################
    minpopx=-1

    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=0.25

#    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
#    mintauv=0.05
#    mintauvneb=0.05
#    maxtauvneberr=0.25
#
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#    
#    mintauv=0.01
#    mintauvneb=0.01
#    maxtauvneberr=999.0
#
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

    minpopx=0.01
    
    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=0.25

#    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
#    mintauv=0.05
#    mintauvneb=0.05
#    maxtauvneberr=0.25
#
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#    
#    mintauv=0.01
#    mintauvneb=0.01
#    maxtauvneberr=999.0
#
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
# 
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

    minpopx=0.05
    
#    mintauv=0.01
#    mintauvneb=0.01
#    maxtauvneberr=0.25
#
##    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
##    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
##    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
##    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
    mintauv=0.05
    mintauvneb=0.05
    maxtauvneberr=0.25
      
    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    OPTS="--rbinstep $rbinstep --v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}_Rs${rbinstep}HLR"
    echo "./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
    fi

#    mintauv=0.05
#    mintauvneb=0.05
#    maxtauvneberr=0.25
#
#    OPTS="--underS06 --v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_underS06_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
#    mintauv=0.01
#    mintauvneb=0.01
#    maxtauvneberr=999.0
#
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
# 
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

    minpopx=0.1
    
#    mintauv=0.01
#    mintauvneb=0.01
#    maxtauvneberr=0.25
#
##    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
##    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
##    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
##    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
    mintauv=0.02
    mintauvneb=0.02
    maxtauvneberr=0.25

    OPTS="--rbinstep $rbinstep --v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}_Rs${rbinstep}HLR"
    echo "./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
    fi

#    mintauv=0.02
#    mintauvneb=0.02
#    maxtauvneberr=0.25
#
#    OPTS="--underS06 --v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_underS06_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#    mintauv=0.01
#    mintauvneb=0.01
#    maxtauvneberr=999.0
#
#    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
# 
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

    minpopx=-1
    mintauv=-1
    mintauvneb=-1
    maxtauvneberr=999.0

    OPTS="--rbinstep $rbinstep --v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    #output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}_Rs${rbinstep}HLR"
    echo "./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $global_opts $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
    fi

#    mintauv=-1
#    mintauvneb=-1
#    maxtauvneberr=999.0
#
#    OPTS="--underS06 --v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    #OPTS="--v_run -1 -R -G --rgbcuts -L $file --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_underS06_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log"
#    if ! $dryrun
#    then
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 4.0 -H ${output}_4HLR.h5 &> ${output}_4HLR.log
#    fi
#
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

done
