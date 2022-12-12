CC = g++
CFLAGS = -Wall -O3 --std=c++11

esegui :
	./es1
compila : es1.cpp es1.o random.o
	$(CC) random.o es1.o -o es1 $(CFLAGS)
es1.o : es1.cpp random.h VectorOperations.h
	$(CC) -c es1.cpp -o es1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
