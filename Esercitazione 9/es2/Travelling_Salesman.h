
#include "random.h"
#include "/Users/alessandro/Libreria/armadillo-11.1.1/include/armadillo"
using namespace std;


int seed[4];
Random rnd;

const int n_city=34;
const int n_config=100;
int in_circle;
double r,x;
double dist=0.,best_dist=0.; //min distance
int index_min;
int i_ran=1,j_ran=1;
int n_steps=1500;

double X[n_city];
double Y[n_city];

int check=0;


arma::Mat<int> A(n_city,n_config);

ofstream Distance;

//pigreco
const double pi=3.1415927;


//functions
void Input(void);
int Cost();
double Norm(double,double,double,double);
void First_Travel();
void Check();
int Min(double* , int );
void Compute_distance(int,int,bool );
void Reorginize();
void Selection_operator();
arma::Mat<int> Crossover(int,int);
bool search_in(int,int*);
arma::Mat<int> Mutation_1(arma::Mat<int>);
arma::Mat<int> Mutation_2(arma::Mat<int>);
void Print_coordinates();









