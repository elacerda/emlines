#!/bin/bash
./SFR.py --minpopx 0.05 --mintauv 0.05 --mintauvneb 0.05 -H SFR_050505.h5 &> SFR_050505.log &
./SFR.py --minpopx 0.20 --mintauv 0.20 --mintauvneb 0.20 -H SFR_202020.h5 &> SFR_202020.log &
./SFR.py --minpopx 0.10 --mintauv 0.10 --mintauvneb 0.10 -H SFR_101010.h5 &> SFR_101010.log 
./SFR.py --minpopx 0.05 --mintauv 0.01 --mintauvneb 0.01 -H SFR_050101.h5 &> SFR_050101.log &
./SFR.py --minpopx 0.10 --mintauv 0.01 --mintauvneb 0.01 -H SFR_100101.h5 &> SFR_100101.log &
./SFR.py --minpopx 0.20 --mintauv 0.01 --mintauvneb 0.01 -H SFR_200101.h5 &> SFR_200101.log 
./SFR.py --spiral --minpopx 0.10 --mintauv 0.10 --mintauvneb 0.10 -H SFR_101010_spiral.h5 &> SFR_101010_spiral.log &
./SFR.py --spiral --minpopx 0.05 --mintauv 0.05 --mintauvneb 0.05 -H SFR_050505_spiral.h5 &> SFR_050505_spiral.log &
./SFR.py --spiral --minpopx 0.05 --mintauv 0.01 --mintauvneb 0.01 -H SFR_050101_spiral.h5 &> SFR_050101_spiral.log  
./SFR.py --spiral --minpopx 0.10 --mintauv 0.01 --mintauvneb 0.01 -H SFR_100101_spiral.h5 &> SFR_100101_spiral.log &
./SFR.py --spiral --minpopx 0.20 --mintauv 0.20 --mintauvneb 0.20 -H SFR_202020_spiral.h5 &> SFR_202020_spiral.log &
./SFR.py --spiral --minpopx 0.20 --mintauv 0.01 --mintauvneb 0.01 -H SFR_200101_spiral.h5 &> SFR_200101_spiral.log 
