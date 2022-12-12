#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

#include "RW.h"



using namespace std;

class Position {

public:

    Position(){
        r_x=0;
        r_y=0;
        r_z=0;
    }
    
    Position(double x, double y, double z){
        r_x=x;
        r_y=y;
        r_z=z;
    }
    
    double Square_Distance(){
        return r_x*r_x+r_y*r_y+r_z*r_z;
    }

    double r_x;
    double r_y;
    double r_z;
    
    protected:
};

int main(){
    
    int n_simulations=N_SIMULATION;
    //int N=1000;
    int n_block=N_BLOCK;
    double direction=0;
    double step=0;
    int n_steps=N_STEPS;
    Position pos;
    double a=1;
    double distance[N_STEPS][N_SIMULATION];
    double block_value[N_STEPS][N_BLOCK];
    
    
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
       
    
    for (int j=0;j<n_simulations;j++){//number of simulation
        pos.r_x=0;
        pos.r_y=0;
        pos.r_z=0;
        for(int i=0;i<n_steps;i++){//RW steps
            direction=rnd.Rannyu(0,3);
            step=rnd.Rannyu(0,2);
            
            if (direction<=1){
                
                if(step<1) pos.r_x=pos.r_x+a;
                if(step>1) pos.r_x=pos.r_x-a;
                
            }
            
            else if (direction>1 and direction<=2){
                
                if(step<1) pos.r_y=pos.r_y+a;
                if(step>1) pos.r_y=pos.r_y-a;
                
            }
            else if (direction>2){
                
                if(step<1) pos.r_z=pos.r_z+a;
                if(step>1) pos.r_z=pos.r_z-a;
                
            }

            distance[i][j]=pos.Square_Distance();
            
        }
        
    }
   
    
    //divide n_simulations in n_block and compute for each block a mean value for each step
    ofstream fout("descrete_RW.txt");
    
    for(int j=0;j<n_steps;j++){ //cycling on steps
        
        for(int i=0;i<n_block;i++){
            
            block_value[j][i]=sqrt(Media_vec(distance,j,i*n_simulations/n_block,(i+1)*n_simulations/n_block));
        }
        fout<<j+1<<" "<<CalcolaMedia(block_value,j)<<" "<<CalcolaErrore(block_value,j)<<endl;
    }

    //----------------------Continue RW-----------------------
    
    double theta=0;
    double phi=0;
    
    for (int j=0;j<n_simulations;j++){//number of simulation
        pos.r_x=0;
        pos.r_y=0;
        pos.r_z=0;
        for(int i=0;i<n_steps;i++){//RW steps
            theta=rnd.Rannyu(0,2*M_PI);
            phi=rnd.Rannyu(0,M_PI);
            
            pos.r_x=pos.r_x+a*sin(phi)*cos(theta);
            pos.r_y=pos.r_y+a*sin(phi)*sin(theta);
            pos.r_z=pos.r_z+a*cos(phi);

            distance[i][j]=pos.Square_Distance();
            
        }
        
    }
    
    ofstream fout2("continue_RW.txt");
    
    for(int j=0;j<n_steps;j++){ //cycling on steps
        
        for(int i=0;i<n_block;i++){
            
            block_value[j][i]=sqrt(Media_vec(distance,j,i*n_simulations/n_block,(i+1)*n_simulations/n_block));
        }
        fout2<<j+1<<" "<<CalcolaMedia(block_value,j)<<" "<<CalcolaErrore(block_value,j)<<endl;
    }
    
    
    fout.close();
    fout2.close();
        return 0;
    }
