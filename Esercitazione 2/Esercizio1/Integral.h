#ifndef __Integral_h__
#define __Integral_h__
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

using namespace std;

double Taylor_exp3(double);
double Taylor_exp2(double);
double Integral( int,int,int,double);
double CalcolaMedia( const vector<double> &, int);
double CalcolaVarianza( const vector<double> & , int );
double CalcolaErrore( const vector<double> & , int );
#endif
