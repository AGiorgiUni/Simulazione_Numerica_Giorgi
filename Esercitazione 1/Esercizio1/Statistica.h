#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "VectorOperations.h"

using namespace std;


 double CalcolaMedia( const vector<double> & v, int n) {

   
     double media=0;
    for(int i=0;i<n;i++){
        
        media=media+v[i];
    }
   
return (double)media/n;

};


 double CalcolaErrore( const vector<double> & v, int n) {

     double accu=0;
    for(int i=0;i<n;i++){
        accu+=(v[i]*v[i]);
    }
     double media=CalcolaMedia(v,n);
return (double)sqrt((accu/n-media*media)/(n-1));

};


double CalcolaSigma( const vector<double> & v, int n) {

    double accu=0;
   for(int i=0;i<n;i++){
       accu+=((v[i]-0.5)*(v[i]-0.5));
   }
    double media=CalcolaMedia(v,n)-0.5;
return (accu/n-media*media);

};

double CalcolaChi( const vector<double> & v, int n, int M) {

    double ni=0;
    double chi=0;
    double min,max=0;
   for(int i=0;i<M;i++){ //Dividing in ranges
       ni=0;

       for(int j=0;j<v.size();j++){ //checking if randomly numbers in vector are in range [i/M,(i+1)/M]
           min=(double)i/M;
           max=(double)(i+1)/M;
           if(min<=v[j] and v[j]<max) {
               ni++;
           }
       }
       cout<<ni<<" numbers found in  ["<<(double)i/M<<","<<(double)(i+1)/M<<"]"<<endl;
       chi=chi+((ni-n/M)*(ni-n/M))/(n/M);
   }
    
    return chi;

};

