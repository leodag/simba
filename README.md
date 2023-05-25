# SIMBA: Similarity-Based Alignment for NMR spectra

This program intends to compare two or more NMR spectra, quantifying how similar they are to one another as a percentage.

## Installation

First, clone this repository:

```
$ git clone https://github.com/leodag/simba.git
```

To be ready to use, you just need to install the dependencies:

```
$ pip install numpy pandas matplotlib
```

## Usage

To use this program, you need to put all the spectra you wish to compare ih tne same directory. Then, you can run it with a command line such as this, for 1H spectra:

```
$ ./simba --directory "data/" --snip-edges -0.3,10.1 --scale-sum --similarity-align 10 \
  --remove-range 7.2,7.32 --snip-edges -0.2,10 --scale-sum --similarity-matrix
```

For 13C spectra, you can use this:

```
$ ./simba --directory "dados/aluno1-13h" --snip-edges -5.1,210.1 --cancel-negative \
--scale-sum --similarity-align 10 --remove-range 76.4,77.6 --snip-edges -5,210 \
--rolling-average 31 --scale-sum --similarity-matrix
```
