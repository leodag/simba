#!/usr/bin/env bash

set -e

run() {
    # $1 is directory and $2 is suffix
    args=(--directory "$1" --snip-edges -0.3,10.1 --scale-sum --similarity-align 10 --remove-range 7.2,7.32 --snip-edges -0.2,10 --scale-sum --similarity-matrix)
    if [[ $# == 3 ]]; then
        args=(--format "$3" "${args[@]}")
    fi
    ./comparator.py "${args[@]}"
    mv out_matrix "results/out_matrix_$2.csv"
    mv out_align "results/out_align_$2.csv"
}

run-1h() {
    args=(--directory "$1" --snip-edges -0.2,11 --scale-sum --similarity-align 10 --similarity-matrix)
    if [[ $# == 3 ]]; then
        args=(--format "$3" "${args[@]}")
    fi
    ./comparator.py "${args[@]}"
    mv out_matrix "results/out_matrix_$2.csv"
    mv out_align "results/out_align_$2.csv"
}

run-1h-no-align() {
    args=(--directory "$1" --snip-edges -0.2,11 --scale-sum --similarity-align 0 --similarity-matrix)
    if [[ $# == 3 ]]; then
        args=(--format "$3" "${args[@]}")
    fi
    ./comparator.py "${args[@]}"
    mv out_matrix "results/out_matrix_$2.csv"
    mv out_align "results/out_align_$2.csv"
}

run-1h-roll-avg() {
    args=(--directory "$1" --snip-edges -0.2,11 --scale-sum --similarity-align 10 --rolling-average 151 --similarity-matrix)
    if [[ $# == 3 ]]; then
        args=(--format "$3" "${args[@]}")
    fi
    ./comparator.py "${args[@]}"
    mv out_matrix "results/out_matrix_$2.csv"
    mv out_align "results/out_align_$2.csv"
}

run-13c() {
    args=(--directory "$1" --snip-edges -5.1,210.1 --scale-sum --similarity-align 10 --remove-range 76.4,77.6 --snip-edges -5,210 --rolling-average 31 --scale-sum --similarity-matrix)
    if [[ $# == 3 ]]; then
        args=(--format "$3" "${args[@]}")
    fi
    ./comparator.py "${args[@]}"
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

# run dados/JanACD janacd
# run dados/Ronald2 ronald2

# run dados/JanTopsp jantopsp TopSpin2
# run dados/JunioTopsp juniotopsp TopSpin2
# run dados/LaisTopsp laistopsp TopSpin2
# run dados/RonaldTopsp ronaldtopsp TopSpin2

# run dados/Gabriel2 gabriel2
# run dados/Ronald3 ronald3
# run dados/Seldon3 seldon3

# run dados/Seldon4Sp seldon4sp TSV

# run dados/LaisTopsp2 laistopsp2 TopSpin2

# run-1h dados/1H-NMR 1h-nmr
# run-1h-roll-avg dados/Vetiver vetiver
# run-1h-roll-avg dados/VetiverIS vetiver-is

# run-13c dados/13C-NMR 13c-nmr

# run-1h dados/David david

# run-1h teste1 teste1 ACD

# run-1h dados/1H-NMR 1h-nmr
run-1h-no-align dados/1H-NMR 1h-nmr-no-align
