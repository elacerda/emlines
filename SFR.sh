#!/bin/bash
#for file in /Users/lacerda/CALIFA/list_SFR/rem_ba_morph/list_0123.txt /Users/lacerda/CALIFA/listOf300GalPrefixes.txt
#dryrun=true
dryrun=false
#for file in /Users/lacerda/CALIFA/listOf298GalPrefixes.txt
for file in /Users/lacerda/CALIFA/listv20_q050.d15a.txt
do
    bname=$(basename $file)
    name=${bname%.txt}
    #minpopx=0.01
    minpopx=0.00
    #mintauv=0.01
    mintauv=0.05
    #mintauvneb=0.01
    mintauvneb=0.05

#    maxtauvneberr=999.0
#    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi
#
#    maxtauvneberr=0.15
#    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi
#
#    ###############
#    ### RGBCuts ###
#    ###############
#    maxtauvneberr=999.0
#    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi
#
#    maxtauvneberr=0.15
#    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi
#
#    #######################
#    ### Residual Filter ###
#    #######################
#    maxtauvneberr=999.0
#    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi
#
#    maxtauvneberr=0.15
#    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi
#
#    ############################
#    ### RGBCuts and Residual ###
#    ############################
#    maxtauvneberr=999.0
#    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
#    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
#    if ! $dryrun
#    then
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
#    time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
#    fi

    #maxtauvneberr=0.15
    maxtauvneberr=0.25
    OPTS="-G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi

    maxtauvneberr=0.5
    OPTS="-G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
done
