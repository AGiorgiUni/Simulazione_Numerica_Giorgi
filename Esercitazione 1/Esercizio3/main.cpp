#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

#include "Statistica.h"

using namespace std;

int main(){
    
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    double d=4; //distance of lines
    double L=3; //length
    //my idea is to generate number beetween 0 and d and cosine of the angle beetween the stick and the horizontal line
    int N_throw=0,N_in=0;
    int passo=100;
    int N=10000;
    int n=N/passo;
    vector<double> P(n);
    vector<double> media_P(n);
    vector<double> errore_P(n);
    double cosine=0.0;
    double x=0.;
    
    for(int j=1;j<=n;j++){
        N_throw=0;
        N_in=0;
        for(int i=0;i<j*n;i++){
            x=d*rnd.Rannyu();
            cosine=rnd.Rannyu();
            N_throw++;
            
            if(x+L*cosine>=d) N_in++;
            
        }
        P[j]=(double)N_throw/N_in;
    }
       
    
    ofstream fout("pi.txt");
    
    for(int i=1;i<=n;i++){
        media_P[i]=CalcolaMedia(P,i+1);
        errore_P[i]=CalcolaErrore(P,i+1);
    }
    for(int i=5;i<=n;i++){
        fout<<(i)*n<<" "<<2*L*media_P[i]/d<<" "<<2*L*errore_P[i]/d<<endl;
        
    }
    
    fout.close();
    
    return 0;
}

