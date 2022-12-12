CC = g++
CFLAGS = -Wall -O3 --std=c++11

esegui :
	./main
compila : main.cpp main.o random.o RW.o
	$(CC) random.o RW.o main.o -o main $(CFLAGS)
main.o : main.cpp random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
RW.o: RW.cpp
	g++ -c RW.cpp -o RW.o
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main seed.out
