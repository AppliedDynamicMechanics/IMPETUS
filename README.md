# **IMPETUS**: **I**nteractive **M**ulti**P**hysics **E**nvironmen**T** for **U**nified **S**imulations

IMPETUS is an object oriented, easy-to-use, high performance, C++ program for three-dimensional simulations of complex physical systems that can benefit a large variety of research areas, especially in cell mechanics. The program implements cross-communication between locally interacting particles and continuum models residing in the same physical space while a network facilitates long-range particle interactions. Message Passing Interface is used for inter-processor communication for all simulations.

A detailed tutorial for IMPETUS is provided in [TUTORIAL.pdf](TUTORIAL.pdf). Documentation and usage examples are available on the [IMPETUS website](http://engr.uconn.edu/~gelyko/impetus.html).

For more information about this project, see “[IMPETUS – Interactive MultiPhysics Environment for Unified Simulations](http://dx.doi.org/10.1016/j.jbiomech.2016.10.042)” (Vi Q. Ha & George Lykotrafits), which was published in the *Journal of Biomechanics* in 2016.

## ***How to run this code***
1. ***Requirements:*** IMPETUS uses OpenMPI for code compilation. Please ensure that OpenMPI is installed and locatable by running the following command in the terminal, which should print the path to your installation of OpenMPI's `mpiCC` compilation module:
    ```
    which mpiCC
    ```

2. ***Compile:*** IMPETUS can be compiled using the provided makefile. Make sure that you're in the `IMPETUS` top directory in your terminal, then run this command:
    ```
    make
    ```

3. ***Run:*** The `simrun` script has been provided as an easy way to execute a compiled IMPETUS program. To invoke this script, you'll need to provide an argument that specifies the number of parallel processes that the program should utilize. This number should not exceed the number of CPU cores in your computer. For example, on a computer with 4 CPUs, you can run:
    ```
    ./simrun 4
    ```
