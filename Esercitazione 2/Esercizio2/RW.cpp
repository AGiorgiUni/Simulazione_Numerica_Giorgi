#include "RW.h"


using namespace std;

double CalcolaMedia( const double M[][N_BLOCK], int n_step) {

  
    double mean=0;
    for(int i=0;i<N_BLOCK;i++){
       
        mean=mean+M[n_step][i];
    }
  
return mean/N_BLOCK;

};

double Media_vec( const double M[][N_SIMULATION], int n_step,int start, int end) {
    double mean=0;
    for(int i=start;i<end;i++){
       
        mean=mean+M[n_step][i];
    }
  
return mean/(end-start);

};

double CalcolaErrore( const double M[][N_BLOCK], int n_step) {

    double accu=0;
   for(int i=0;i<N_BLOCK;i++){
       accu+=(M[n_step][i])*(M[n_step][i]);
   }
    double media=CalcolaMedia(M,n_step);
return sqrt((accu/N_BLOCK-media*media)/(N_BLOCK-1));

};

