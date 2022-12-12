CC = g++
CFLAGS = -Wall -O3 --std=c++11 

esegui:
	./Travelling_Salesman.exe
compila : Travelling_Salesman.cpp Travelling_Salesman.o random.o
	$(CC) random.o Travelling_Salesman.o -o Travelling_Salesman.exe $(CFLAGS)
Travelling_Salesman.o : Travelling_Salesman.cpp random.o Travelling_Salesman.h 
	$(CC) -c Travelling_Salesman.cpp -o Travelling_Salesman.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o Travelling_Salesman.exe 
