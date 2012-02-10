#!/bin/bash
gcc sar_simulator.c -lfftw3 -lm -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -march=native -g -o sar_simulator
