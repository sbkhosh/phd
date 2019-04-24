#!/bin/bash
#$ -cwd -V
#$ -l h_rt=00:15:00
#$ -l np=16
mpirun inc3d_polaris
