#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "random.h"

#include "Economy.h"

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
    
   
    double S0=100;
    double t=0;
    double T=1;
    double K=100;
    double r=0.1;
    double sigma=0.25;
    int N=10000;
    double gauss=0;
    int M =100;
    int steps=N/M;
    //vector<double> x(10000);
    double x;
    vector<double> C(N);
    vector<double> P(N);
    vector<double> mediaC(steps);
    vector<double> erroreC(steps);
    vector<double> mediaP(steps);
    vector<double> erroreP(steps);
   
    //calculating S(t) directly using the formula
    for (int i=0;i<N;i++){
        gauss=rnd.Gauss(0,T);
        x=S_t(T,r,sigma,S0,gauss);
        C[i]=exp(-r*T)*max(0,x-K);
        P[i]=exp(-r*T)*max(0,K-x);
    } 
    
    ofstream fout("dati.txt");
    
    for (int i=0;i<steps;i++){
        mediaC[i]=CalcolaMedia(C,(i+1)*M);
        mediaP[i]=CalcolaMedia(P,(i+1)*M);
        erroreC[i]=CalcolaErrore(C,(i+1)*M);
        erroreP[i]=CalcolaErrore(P,(i+1)*M);
        
        fout<<(i+1)*M<<" "<<mediaC[i]<<" "<<erroreC[i]<<" "<<mediaP[i]<<" "<<erroreP[i]<<endl;
    }
    
    fout.close();
      //reset vectors of put and call
    fill(C.begin(), C.end(), 0);
    fill(P.begin(), P.end(), 0);
    fill(mediaC.begin(), mediaC.end(), 0);
    fill(mediaP.begin(), mediaP.end(), 0);
    fill(erroreC.begin(), erroreC.end(), 0);
    fill(erroreP.begin(), erroreP.end(), 0);
    
    double step=(double)1/100;
    
    double S_i=S0;
    
    double tf=t+step;
    
    for (int i=0;i<N;i++){
        t=0;
        tf=t+step;
        S_i=S0;
        while(tf<T){
            gauss=rnd.Gauss(0,1);
            S_i=S_t_BySteps(tf,t,r,sigma,S_i,gauss);
            t=t+step;
            tf=tf+step;
            
        }
        x=S_i;
        
        C[i]=exp(-r*T)*max(0,x-K);
        P[i]=exp(-r*T)*max(0,K-x);
    }
    
    ofstream fout2("dati2.txt");
    
    for (int i=0;i<steps;i++){
        mediaC[i]=CalcolaMedia(C,(i+1)*M);
        mediaP[i]=CalcolaMedia(P,(i+1)*M);
        erroreC[i]=CalcolaErrore(C,(i+1)*M);
        erroreP[i]=CalcolaErrore(P,(i+1)*M);
        
        fout2<<(i+1)*M<<" "<<mediaC[i]<<" "<<erroreC[i]<<" "<<mediaP[i]<<" "<<erroreP[i]<<endl;
    }
    
    fout2.close();
    
    return 0;
    }

