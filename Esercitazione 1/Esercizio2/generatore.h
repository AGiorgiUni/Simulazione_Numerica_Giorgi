#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "VectorOperations.h"

using namespace std;


 double ExpNumber( double lambda, double x) {

     if(0<=x<=1.){
         return 1-exp(-lambda*x);
     }
     else {cout<<"Error"<<endl;
         return -1;
     }

};


double LorentzNumber( double mu, double x, double gamma) {

    if(0<=x<=1.){
        return atan((x-mu)/gamma)/M_PI+0.5;
    }
    else {cout<<"Error"<<x<<endl;
        return -1;
    }

};
