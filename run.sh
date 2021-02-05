#!/usr/bin/env bash

run() {
    # $1 is directory and $2 is suffix
    ./comparator.py --directory "$1" --snip-edges -0.3,10.1 --scale-sum --similarity-align 10 --remove-range 7.2,7.32 --snip-edges -0.2,10 --scale-sum --similarity-matrix
    mv out_matrix "results/out_matrix_$2.csv"
    mv out_align "results/out_align_$2.csv"
}

# run dados/Jan jan
# run dados/Gabriel gabriel
# run dados/Junio junio
# run dados/Lais lais
# run dados/Lais2 lais2
# run dados/Ronald ronald
# run dados/Seldon seldon

run dados/JanACD janacd
run dados/Ronald2 ronald2
# run dados/JanTopsp jantopsp
