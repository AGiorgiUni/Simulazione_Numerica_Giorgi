#include "Economy.h"
#include "random.h"
#include <vector>

using namespace std;

double W(double x,double tf){
    
    if(tf==0){
        return 0;
    }
    
    else{
        
        return x*sqrt(tf);
    }
    
}


double S_t(double t,double mu,double sigma,double S0,double x){
    
    return S0*exp((mu-0.5*sigma*sigma)*t+sigma*x);
    
}

double S_t_BySteps(double tf,double ti,double mu,double sigma,double S_i,double x){
    
    return S_i*exp((mu-0.5*sigma*sigma)*(tf-ti)+sigma*x*sqrt(tf-ti));
    
}

double N(Random rnd){
    
    return 0.5*(1+erf(rnd.Gauss(0,1)/(sqrt(2))));
}

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
       accu=accu+(v[i]*v[i]);
   }
    double media=CalcolaMedia(v,n);
return sqrt((accu/n-media*media)/(n-1));

};


double max(double a, double b){
    if(a>=b) return a;
    else return b;
}

