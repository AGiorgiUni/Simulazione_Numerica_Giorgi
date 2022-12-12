#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "function.h"

using namespace std;

int main()
{
    
    
    Input();
    Ene.open("output.ene.0", ios::app);
    mu=0.5;
    sigma=0.5;
    cout<<endl<<"Mu range = ["<<mu<<","<<mu+step_4_temp*0.05<<"]"<<endl;
    cout<<endl<<"Sigma range = ["<<sigma<<","<<mu+step_4_temp*0.05<<"]"<<endl;
    while(temp>0.001){
        //change temp means changing probability of metropolis
        //mu and sigma can't be random or problem in evaluating energy
        for (int i=0;i<step_4_temp;i++){
            //sigma_try=0.55;
            mu_try = 0.5+i*0.05;
            for (int j=0;j<step_4_temp;j++){
                sigma_try = 0.5+j*0.05;
                H = Hamiltonian(false, mu, sigma);
                H_try=Hamiltonian(false, mu_try, sigma_try);
                p=minimo(1,exp(-beta*(H_try-H)));
                attempted++;
                if(rnd.Rannyu()<=p){
                    mu=mu_try;
                    sigma=sigma_try;
                    accepted++;
                }
        
            }

        }
        
        cout<<"Temperature = "<< temp<<endl;
        cout << "Acceptance rate " << (double)accepted/attempted << endl << endl;
        cout << "----------------------------" << endl << endl;
        Print_1(temp, mu, sigma);
        
        temp=temp/10;
        beta=1/temp;
    }
    
    H = Hamiltonian(true, mu, sigma);
    Histo(mu, sigma);

    Ene.close();
  return 0;
}

void Input(void)
{
    
    ifstream Primes, Seed;
    ifstream open("dati.dat");
    //Read seed for random numbers
      int p1, p2;
      Primes.open("Primes");
      Primes >> p1 >> p2 ;
      Primes.close();
    
        open >> restart;
        //open>>delta;

      if(restart) Seed.open("seed.out");
      else Seed.open("seed.in");
      Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
      rnd.SetRandom(seed,p1,p2);
      Seed.close();
 
    
    cout << "The program uses h_cut=1 and m=1 " << endl;

    open >> temp;
    cout<<"Starting temp = "<<temp<<endl;
    beta=1/temp;
    open.close();

  return;
}


double Energy(double x,double mu, double sigma){
    
    return -0.5*Derivata2(x,mu,sigma)+Potential(x);
    
}
double Potential(double x){
    
    return pow(x,4)-5*pow(x,2)/2;
    
}


double Derivata2(double x,double mu, double sigma){
    return (pow((x-mu),6)/(4*pow(sigma,8)))*exp(-pow(x-mu,2)/(2*sigma*sigma))+(pow((x+mu),6)/(4*pow(sigma,8)))*exp(-pow(x+mu,2)/(2*sigma*sigma));
}

double Error( int n, double sum, double sum2) {
    if (n == 1) return 0.;
    else return sqrt((sum2 / double(n) - pow(sum / double(n), 2)) /double(n - 1));
}

double minimo(double a, double b){
    
    if(a<=b) return a;
    else return b;
}

void Print_1(double temp, double mu, double sigma){
    param.open("output.param.0",ios::app);
    H=Hamiltonian(false,mu,sigma);
    param << temp<<" "<<mu<<" "<<sigma<<" " << H << " "<< err_H << endl;
    param.close();
}

double Hamiltonian(bool print, double mu, double sigma){
    //calculate energy value for mu and sigma
    
    p=0.0;
    sum=0.;
    sum2=0.;
    x=rnd.Gauss(mu,sigma);
    for(int i=0;i<n_blocks;i++){
        block_avg=0;
        for (int j=0;j<n_steps;j++){
            delta=mu+rnd.Rannyu(-2,2)*sigma; //metodo per passare da una gaussiana all'altra
            if(rnd.Rannyu()<0.5) delta=(-1)*delta;
            y=x-delta;
            //attempted++;
            p_old=Eval(x,mu,sigma); //probability density
            p_new=Eval(y,mu,sigma);
            p=minimo(1,p_new/p_old);
            if(rnd.Rannyu()<=p){
                x=y;
                //accepted++;
            }
            block_avg+=Energy(x,mu,sigma);
        }
        
        block_avg /= n_steps;
        sum += block_avg;
        sum2 +=pow(block_avg, 2);
        
        stima_H=sum / double(i + 1);
        err_H=Error(i + 1, sum, sum2);
        
        if(print){
            Ene << i + 1 << " " << stima_H << " "<< err_H << endl;

        }
        
        /*cout << "Block number " << i+1 << endl;
        //Averages(i,blk_avg[i],blk_avg2[i]);
        cout << "Acceptance rate " << (double)accepted/attempted << endl << endl;
        cout << "----------------------------" << endl << endl;*/
        
      }
    
    return stima_H;
}

void Histo(double mu, double sigma){
    
    double dr=(double)6/100; //sampling in [-2,2]
    x=rnd.Gauss(mu,sigma);
    for (int j=0;j<n_steps;j++){
        delta=mu+rnd.Rannyu(-2,2)*sigma; //metodo per passare da una gaussiana all'altra
        if(rnd.Rannyu()<0.5) delta=(-1)*delta;
        //y=x-delta;
        y=rnd.Rannyu(mu-3*sigma,mu+3*sigma);
        //attempted++;
        p_old=Eval(x,mu,sigma); //probability density
        p_new=Eval(y,mu,sigma);
        p=minimo(1,p_new/p_old);
        if(rnd.Rannyu()<=p){
            x=y;
            //accepted++;
        }
        bin_index = (int)(3/dr+(x/dr)); // posizione nell'istogramma
        histo[bin_index ] += 1;

    }
    
        ofstream Histo;
        Histo.open("output.Histo.0",ios::app);

        for (int i=0;i<n_bins;i++){
            Histo << -2+i*dr<<" "<<(double)histo[i]/n_steps<< endl;
        }
        Histo.close();
    
    
    return ;
}


