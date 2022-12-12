#ifndef __RW_h__
#define __RW_h__
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

using namespace std;

const int N_BLOCK=100;
const int N_SIMULATION=10000;
const int N_STEPS=100;

double CalcolaMedia( const double[][N_BLOCK], int);
double Media_vec( const double [][N_SIMULATION], int ,int , int );
double CalcolaErrore( const double [][N_BLOCK], int );
#endif
