# acse-6-individual-assignment-acse-efb119
acse-6-individual-assignment-acse-efb119 created by GitHub Classroom

## Introduction
This project required the implementation of MPI into Conway's Game of Life in C++, the program outputs 2 data files per core. Both are used for postprocessing done in the animation.py file.

## Installation:
Clone the repository from Github by either:
* Using command line:
`https://github.com/acse-2019/acse-6-individual-assignment-acse-efb119.git`
* Downloading the repository as a .zip

## Usage:
Using Unix command line arguments for the compiler from within the library directory, run:

```mpicxx -o life LifeLocal.cpp```

Then run:

```mpiexec -n 4 ./life```

to execute the program. (where the ```-n 4``` is the number of cores in this case 4 to be used).

For the postprocessing run:

```python animation.py```

Some issues with this postproccesing program it works mainly for 4 cores and likes having processor grids to have the same number of columns. If time permitted would fix so can work for any shape and grid sizing!

## Input and Output
The LifeLocal.cpp file has 4 variables in the code that can be changed to change output, Rows, Cols, Period and Periodic. The file will output a file that has the processor dimensions and a file that has the game of life grid for each processor.

The animation.py file will ask for user input in terminal about how many cores were used to run the .cpp file. This file saves a gif of game of life.
