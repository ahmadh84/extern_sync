#!/bin/bash

MYOLDPWD=${PWD};
for f in {cars2,cars3,cars6,cars7,cars8,cars9,marple1,marple3,marple5,marple8,marple10,marple11,marple13,bear01,bear02,cats02,cats04,cats05,cats07,ducks01,horses01,horses03,horses06,lion02,meerkats01,people04,people05,rabbits01,rabbits05}; do 
    cd TrainingSet/Data/${f}/;
    mogrify -format ppm *.jpg;
    cd ${MYOLDPWD};
done;

