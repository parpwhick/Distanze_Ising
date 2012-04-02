COPTS = -Wall -fopenmp
COPTS += -DRIDUZIONE
#COPTS += -g 
COPTS += -O3

files=*.cpp *.h Makefile2

ising: ising_disordinato.o distanze.o partizioni.o init_functions.o rand55.o
	g++ -o ising ising_disordinato.o distanze.o partizioni.o init_functions.o rand55.o -lgomp

rand55.o: rand55.cpp rand55.h
	g++ ${COPTS} -c rand55.cpp

ising_disordinato.o: ising_disordinato.cpp strutture.h
	g++ ${COPTS} -c ising_disordinato.cpp

main.o: main.cpp strutture.h 
	g++ ${COPTS} -c main.cpp

translation.o: strutture.h translation.cpp
	g++ ${COPTS} -c translation.cpp

init_functions.o: strutture.h init_functions.cpp
	g++ ${COPTS} -c init_functions.cpp

distanze.o: strutture.h distanze.cpp
	g++ ${COPTS} -c distanze.cpp

partizioni.o: strutture.h partizioni.cpp
	g++ ${COPTS} -c partizioni.cpp

clean:
	rm -f *.o
