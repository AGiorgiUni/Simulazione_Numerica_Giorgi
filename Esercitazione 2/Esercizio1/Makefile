CC = g++
CFLAGS = -Wall -O3 --std=c++11 

esegui :
	./main
compila : main.cpp main.o random.o Integral.o
	$(CC) random.o Integral.o main.o -o main $(CFLAGS)
main.o : main.cpp random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
Integral.o: Integral.cpp
	g++ -c Integral.cpp -o Integral.o
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main seed.out
