
//Random numbers
#include "random.h"

int seed[4];
Random rnd;
int restart;

double x=0,mu=0,sigma=0,delta=0.;
double y=0.;
double temp=0.,beta=0.;
int attempted=0, accepted=0;

const int n_steps=15000;
const int n_blocks=20;
const int step_4_temp=20;

double p_old,p_new,p;

double block_avg;

double sum=0.;
double sum2=0.;

double stima_H, err_H;

double H=0.,H_try=0.;
double mu_try=0.,sigma_try=0.;

const int n_bins=101;
int histo[n_bins];
int bin_index=0;

std::ofstream Ene;
std::ofstream param;

//functions

double Eval(double x,double mu, double sigma){
    return pow((exp(-pow(x-mu,2)/(2*sigma*sigma))+exp(-pow(x+mu,2)/(2*sigma*sigma))),2);
}
double Generator(double, double, double);
void Averages(int iblk);
void Input();
double Energy(double ,double, double);
double Potential(double);
double Derivata2(double,double, double);
void Averages(int,double,double);
double Error(int,double, double);
double minimo(double, double);
double Hamiltonian(bool, double, double);
void Print_1(double, double, double);
void Histo(double , double );


