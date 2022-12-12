#include "Integral.h"
using namespace std;

double Taylor_exp3(double x){
    
    return abs(M_PI/2*(1-pow(M_PI*x/2,2)/8+pow(M_PI*x/2,4)/24));
}

double Taylor_exp2(double x){
    
    return abs(M_PI/2*(1-pow(M_PI*x/2,2)/8));
}


 double Integral( int b, int a,int n,double x) {

     double f=0;
    
     int passi=0;
     while(passi<n){
         
         f=f+(M_PI/2)*cos(M_PI*x/2);
         
         passi++;
     }
     return (b-a)*f/n;
};


double CalcolaMedia( const vector<double> & v, int n) {

  
    double media=0;
   for(int i=0;i<n;i++){
       
       media=media+v[i];
   }
  
return (double)media/n;

};

double CalcolaVarianza( const vector<double> & v, int n) {

    double accu=0;
   for(int i=0;i<n;i++){
       accu+=((v[i])*(v[i]));
   }
    double media=CalcolaMedia(v,n);
return (accu/n-media*media);

};


 double CalcolaErrore( const vector<double> & v, int n) {

     double accu=0;
    for(int i=0;i<n;i++){
        accu+=(v[i]*v[i]);
    }
     double media= CalcolaMedia(v,n);
return (double)sqrt((accu/n-media*media)/(n-1));

};

