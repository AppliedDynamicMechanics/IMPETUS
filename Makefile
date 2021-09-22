# Change the path to openmpi:
MPI_PATH = /home3/shared_bin/openmpi/bin/

# Change the path to the simulation engine:
IMPETUS_DIR=lib

CC = $(MPI_PATH)mpiCC

bin/vamde-run.exe: object/main.o
	$(CC) object/main.o  -I$(IMPETUS_DIR)/include -o bin/impetus-run.exe -L$(IMPETUS_DIR)/ -lvamde-vamde-p
	
object/main.o: src/main.cpp
	$(CC) -c -w  -I$(IMPETUS_DIR)/include $< -o $@ 

clean:
	rm -f bin/* object/*
