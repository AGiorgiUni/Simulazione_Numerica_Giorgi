CC = g++
CFLAGS = -Wall -O3 --std=c++11 

esegui:
	mpiexec -np 10 ./test file.txt
compila : Travelling_Salesman.cpp random.o
	mpicxx Travelling_Salesman.cpp random.o -o test
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o test 
