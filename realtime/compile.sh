#!/bin/bash
gcc sar_simulator.c -lfftw3 -lm -O3 -march=native -g -o sar_simulator
