CC = g++
CFLAGS = -Wall -O3 --std=c++11

esegui:
	./main.exe
compila : main.cpp main.o random.o
	$(CC) random.o main.o -o main.exe $(CFLAGS)
main.o : main.cpp random.h function.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe 
