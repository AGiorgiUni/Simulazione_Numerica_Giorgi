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
    vector<double> x(1);
    int passo=1000;
    int N=100000;
    int passi=N /passo;
    vector<double> media(passi);
    vector<double> errore(passi);
    vector<double> varianza(passi);
    
    
    for (int i=0;i<N ;i++){
        x.push_back(rnd.Rannyu());
    }
    
    
   
    for(int i=1;i<passi+1;i++){
        
        media[i-1]=CalcolaMedia(x,i*passo);
        errore[i-1]=CalcolaErrore(x,i*passo);
        
    }
    
    ofstream fout("dati.txt");
    
    for(int i=0;i<passi;i++){
        fout<<(i+1)*passo<<" "<<media[i]<<" "<<errore[i]<<endl;
        
    }
    
    fout.close();
    
    
for (int i=0;i<=x.size() ;i++){
          x[i]=(x[i]-0.5)*(x[i]-0.5);
        
    }
    
    for(int i=1;i<passi+1;i++){
        
        media[i-1]=CalcolaMedia(x,i*passo)-(double)1/12;
       errore[i-1]=CalcolaErrore(x,i*passo);
    
    }
     
    ofstream fout2("dati2.txt");
    
    for(int i=0;i<passi;i++){
        fout2<<(i+1)*passo<<" "<<media[i]<<" "<<errore[i]<<endl;
        
    }
    
    fout2.close();
    
    int M=100;
    int n=10000;
    vector<double> x2(n);
    vector<double> chi(M);
   
    
    
    for(int i=0;i<100;i++){
        for (int j=0;j<x2.size() ;j++){ //loading new data
            x2[j]=rnd.Rannyu();
        }
        chi[i]=CalcolaChi(x2,n,M);
        
    }
    
    ofstream fout3("Dati3.txt");
    
    for(int i=0;i<chi.size();i++){
       
        fout3<<i+1<<" "<<chi[i]<<endl;
        
    }
    
    fout3.close();
    
    return 0;
}
