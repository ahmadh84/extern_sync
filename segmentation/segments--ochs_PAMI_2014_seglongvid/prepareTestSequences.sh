#!/bin/bash

MYOLDPWD=${PWD};
for f in {cars1,cars4,cars5,cars10,marple2,marple4,marple6,marple7,marple9,marple12,people1,people2,tennis,camel01,cats01,cats03,cats06,dogs01,dogs02,farm01,goats01,horses02,horses04,horses05,lion01,people03,giraffes01,rabbits02,rabbits03,rabbits04}; do
    cd TestSet/Data/${f}/;
    mogrify -format ppm *.jpg;
    cd ${MYOLDPWD};
done;

