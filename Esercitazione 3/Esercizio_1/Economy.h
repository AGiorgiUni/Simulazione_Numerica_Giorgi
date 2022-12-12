#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

using namespace std;

double W(double ,double ,double );

double S_t(double ,double ,double ,double ,double );
double N(Random );
double CalcolaMedia( const vector<double> & , int );
double CalcolaErrore( const vector<double> & , int );
double CalcolaVarianza( const vector<double> & , int );
double max(double,double);
double S_t_BySteps(double ,double ,double ,double ,double ,double );
