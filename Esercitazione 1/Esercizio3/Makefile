CC = g++
CFLAGS = -Wall -O3 --std=c++11

esegui :
	./main
compila : main.cpp main.o random.o
	$(CC) random.o main.o -o main $(CFLAGS)
main.o : main.cpp random.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main seed.out
