COPTS= -O3 -march=native 
# debugging symbols
#COPTS+=-g2
# debugging the STL library
#COPTS+=-D_GLIBCXX_DEBUG
# profiling
#COPTS+= -pg
# STL profiling
#COPTS+= -D_GLIBCXX_PROFILE
COPTS +=-Wall -fopenmp
COPTS +=-std=c++0x
LINK_OPTS = -lm -lgomp
files=*.cpp *.h Makefile Doxyfile
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

files=*.cpp *.h Makefile

ising: ising_disordinato.o init_functions.o distance.o partizioni.o adj_handler.o rand_mersenne.o sequence_partitions.o
	g++ -o ising ising_disordinato.o distance.o partizioni.o init_functions.o adj_handler.o rand_mersenne.o sequence_partitions.o -lgomp

ising_disordinato.o: ising_disordinato.cpp strutture.h distance.h partizioni.h
	g++ ${COPTS} -c ising_disordinato.cpp

main.o: main.cpp strutture.h 
	g++ ${COPTS} -c main.cpp

translation.o: strutture.h translation.cpp
	g++ ${COPTS} -c translation.cpp

init_functions.o: strutture.h init_functions.cpp
	g++ ${COPTS} -c init_functions.cpp

distance.o: strutture.h distance.cpp distance.h
	g++ ${COPTS} -c distance.cpp

partizioni.o: strutture.h partizioni.cpp adj_handler.h partizioni.h
	g++ ${COPTS} -c partizioni.cpp

adj_handler.o: adj_handler.cpp adj_handler.h
	g++ ${COPTS} -c adj_handler.cpp
	
rand_mersenne.o: rand_mersenne.cpp rand_mersenne.h
	g++ ${COPTS} -c rand_mersenne.cpp

sequence_partitions.o: sequence_partitions.cpp
	g++ ${COPTS} -c sequence_partitions.cpp

clean: clean_temp_files
	rm -f *.o

clean_temp_files:
	rm -vf *.bin
	rm -vf *.txt

zip: ${files} ${file_supporto}
	zip -9 prog_distanze.zip ${files} ${file_supporto}

arch: ${files}
	mkdir distanze_ising
	cp ${files} distanze_ising
	tar -cvzf prog_distanze_ising.tar.gz distanze_ising
	rm -fr distanze_ising
