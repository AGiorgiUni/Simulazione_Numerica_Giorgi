#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

#include "Integral.h"

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
       
    int M=10000;
    int N=1000;
    int n_block=100;
    
    vector<double> Integral_mean(n_block);
    vector<double> Error(n_block);
    double y=0;
    vector<double> x(M);
   // double interval=(double)1/M;

    ofstream fout("dati_unif.txt");

    for(int k=0;k<M;k++){ //loading x M times with Integral values
        y=rnd.Rannyu();
        x[k]=Integral(1,0,N,y);
    }
   
    for(int j=0;j<n_block;j++){
        //(j+1)*times=(j+1)*M/n_block Intergal values for every block
        Integral_mean[j]=CalcolaMedia(x,(j+1)*M/n_block);
        Error[j]=CalcolaErrore(x,(j+1)*M/n_block);
        fout<<j+1<<" "<<Integral_mean[j]<<" "<<Error[j]<<endl;
    }
   
    
    ofstream fout2("dati_taylor3.txt"); //fourth order taylor exp
    double f_ip=0;
    double y_min=0;
    double y_max=Taylor_exp3(0);
    for(int k=0;k<M;k++){ //loading x M times with Integral values
        y=rnd.Rannyu();
        
        f_ip=rnd.Rannyu(y_min,y_max); //generating number in [0,max] of Taylor expansion function
        while(f_ip>Taylor_exp3(y)){ //accept/reject
            y=rnd.Rannyu();
        }
        x[k]=Integral(1,0,N,y);
    }
    
    for(int j=0;j<n_block;j++){
        Integral_mean[j]=CalcolaMedia(x,(j+1)*M/n_block);
        Error[j]=CalcolaErrore(x,(j+1)*M/n_block);
        fout2<<j+1<<" "<<Integral_mean[j]<<" "<<Error[j]<<endl;
    }
    
    ofstream fout3("dati_taylor2.txt"); //second order taylor exp
    f_ip=0;
    y_min=0;
    y_max=Taylor_exp2(1);
    double r=0;
    for(int k=0;k<M;k++){
        y=rnd.Rannyu();
        
        r=rnd.Rannyu();
        while(r>Taylor_exp2(y)/y_max){
            y=rnd.Rannyu();
        }
        x[k]=Integral(1,0,N,y);
    }
    
    for(int j=0;j<n_block;j++){
        //(j+1)*times=(j+1)*M/n_block Intergal values for every block
        Integral_mean[j]=CalcolaMedia(x,(j+1)*M/n_block);
        Error[j]=CalcolaErrore(x,(j+1)*M/n_block);
        fout3<<j+1<<" "<<Integral_mean[j]<<" "<<Error[j]<<endl;
    }
   
   
    
    fout.close();
    fout2.close();
    fout3.close();
   
        return 0;
    }
