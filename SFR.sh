#!/bin/bash
#dryrun=true
dryrun=false
#for file in /Users/lacerda/CALIFA/listOf298GalPrefixes.txt
#for file in debug.txt
for file in /Users/lacerda/CALIFA/listv20_q050.d15a.txt
do
    bname=$(basename $file)
    name=${bname%.txt}
    
    minpopx=-1
    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=0.25

#    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    
    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi

    minpopx=-1
    mintauv=0.05
    mintauvneb=0.05
    maxtauvneberr=0.25

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
    
    minpopx=-1
    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=999.0

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi

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
    
    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi

    minpopx=0.01
    mintauv=0.05
    mintauvneb=0.05
    maxtauvneberr=0.25

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
    
    minpopx=0.01
    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=999.0

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
 
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

    minpopx=0.05
    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=0.25

#    OPTS="-L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="--rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
#    OPTS="-R --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    
    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi

    minpopx=0.05
    mintauv=0.05
    mintauvneb=0.05
    maxtauvneberr=0.25

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
    
    minpopx=0.05
    mintauv=0.01
    mintauvneb=0.01
    maxtauvneberr=999.0

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
 
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

    minpopx=-1
    mintauv=-1
    mintauvneb=-1
    maxtauvneberr=999.0

    OPTS="--v_run -1 -R -G --rgbcuts -L $file --spiral --minpopx $minpopx --mintauv $mintauv --mintauvneb $mintauvneb"
    output="SFR_${minpopx}_${mintauv}_${mintauvneb}_${maxtauvneberr}_resid_rgbcuts_${name}"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr -H ${output}.h5 &> ${output}.log &"
    echo "./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log"
    if ! $dryrun
    then
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr               -H ${output}.h5      &> ${output}.log &
        time ./SFR.py $OPTS --maxtauvneberr $maxtauvneberr --rbinfin 5.0 -H ${output}_5HLR.h5 &> ${output}_5HLR.log
    fi
 
    #################################################################################################################
    #################################################################################################################
    #################################################################################################################

done
