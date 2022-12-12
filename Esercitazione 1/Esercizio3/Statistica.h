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
return (double)sqrt(accu/n-media*media)/(n-1);

};

double CalcolaVarianza( const vector<double> & v, int n) {

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
   for(int i=0;i<M;i++){ //condizione sulla divisione di [0,M]
       ni=0;
       //cout<<i<<"   "<<ni<<endl<<endl;
       for(int j=0;j<v.size();j++){ //faccio scorrere il vettore e controllo se x è interno
           //cout<<" Trovato tra  "<<(double)i/M<<" e "<<(double)(i+1)/M<<endl;
           if((double)i/M<v[j]<(double)(i+1)/M) {
               //cout<<" Trovato tra  "<<(double)i/M<<" e "<<(double)(i+1)/M<<endl;
               ni++;
           }
       }
       cout<<i<<"   "<<ni<<endl;
       chi=chi+((ni-n/M)*(ni-n/M))/(n/M);
   }
    
    return chi;

};

