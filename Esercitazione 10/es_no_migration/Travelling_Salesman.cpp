#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "mpi.h"
#include "random.h"
#include <armadillo>
#include <string.h>



using namespace std;

bool MIGRATION=false;

int seed[4];
Random rnd;

const int n_city=50;
const int n_config=100;
int in_circle;
double r,x;
double dist=0.,best_dist=0.; //min distance
int index_min;
int i_ran=1,j_ran=1;
int n_steps=10000;

string states[n_city];
string capitals[n_city];
double X[n_city];
double Y[n_city];

int check=0;

ofstream coordinates;

arma::Mat<int> A(n_city,n_config);

ofstream Distance;

//pigreco
const double pi=3.1415927;


//functions
void Input(int);
int Cost();
double Norm(double,double,double,double);
void First_Travel();
void Check();
int Min(double* , int );
void Compute_distance(int,int,bool );
void Reorganize_vec(int*);
void Reorganize();
void Selection_operator();
arma::Mat<int> Crossover(int,int);
bool search_in(int,int*);
arma::Mat<int> Mutation_1(arma::Mat<int>);
arma::Mat<int> Mutation_2(arma::Mat<int>);
void Print_coordinates();
void Print_coordinates_city(int);

int main(int argc, char* argv[])
{
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    
    int h = 0;
    int giver = 0;
    int receiver = 1;
    int msg[n_city];
    int itag=1;
    
    if(rank==0) Distance.open("output.distance.dat",ios::app);
    
    Input(rank);//Inizialization
    First_Travel();//creating travels
    Check();
    double tstart = MPI_Wtime();
    for(int step=0;step<n_steps;step++){
        if(MIGRATION){ //Concept: a rank decides who gives and who receives and then giver broadcast its best sequence
            if(step%20==0){ //migration after 20 steps
                if(rank == 0){
                    giver = (int)rnd.Rannyu(0,size);
                    receiver = (int)rnd.Rannyu(0,size);
                    while(receiver==giver){
                        receiver = (int)rnd.Rannyu(0,size);
                    }
                    cout << endl << "Migration from " << giver << " to " << receiver << endl;
                    
                }
                
                MPI_Bcast(&giver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
                MPI_Bcast(&receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
                
                if(rank == giver){
                    for(int i = 0; i< n_city; i++){
                        msg[i] = A(i,index_min);
                    }
                    MPI_Send(msg,n_city,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD);
                    
                }
                
                if(rank == receiver){
                    MPI_Recv(msg,n_city,MPI_INTEGER,giver,itag,MPI_COMM_WORLD,&stat);
                    Reorganize_vec(msg);
                    Check();
                }
            }
        }
        //same code as the others
        if(rnd.Rannyu()<0.1) A=Mutation_1(A);
        if(rnd.Rannyu()<0.15) A=Mutation_2(A);
        Check();//check boundary conditions
        index_min=Cost();
        if(rank==0) Compute_distance(step,index_min, true);
        else Compute_distance(step,index_min, false);
        Reorganize();
        Check();
        for(int j=1;j<n_config-1;j++){
            if(rnd.Rannyu()<0.5){
                Selection_operator();
                A=Crossover(i_ran,j);
            }
        }
        Check();
    }
    
    
    index_min=Cost();
    if(rank==0) Compute_distance(n_steps,index_min, true);
    else Compute_distance(n_steps,index_min, false);
    //A.print();
    double tend = MPI_Wtime();
    double dt = tend - tstart;
    cout<<" Time with "<<size<< " process = "<<dt<<endl;
    if (rank==0) Print_coordinates();
    Distance.close();
    coordinates.close();
    MPI_Finalize();
  return 0;
}


void Input(int rank)
{
    ifstream ReadInput, ReadConf, ReadVelocity, Seed;
    
    if(rank==0){
        cout << "The Traveling Salesman Problem       " << endl<<endl;
        
        cout << "Using L1 norm: L1(x)=|x(i+1)-x(i)|" << endl;
    }

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
    
    ReadInput.open("input.in");

    ReadInput >> in_circle;
    if(rank==0){
        if(in_circle==1){
            cout << "Cities are located on a circumference" <<endl;
        }
        if(in_circle==0) {
            cout << "Cities are located in a square" << endl;
        }
        if(in_circle==3) {
            cout << "Cities are American capitals" << endl;
        }
    }
  
   if(rank==0) cout << "Number of cities = " << n_city << endl;

    ReadInput >> r;
    if(in_circle==1){
        cout << "Circumference radius = " <<r<<endl;
    }
    if(in_circle==0) {
        cout << "Side of square = "<<r << endl;
    }

    

    if(rank==0) cout << "Number of possible walks = " << n_config << endl;
    
  ReadInput.close();
    
    if(in_circle==1){
        for(int i=0;i<n_city;i++){ //filling the first culumn
            X[i]=r*rnd.Rannyu(-1,1);
            Y[i]=sqrt(r*r-X[i]*X[i]);
            if(rnd.Rannyu()<0.5) Y[i]=Y[i]*(-1);
            check+=i;
        }
    }
    if(in_circle==0){
        for(int i=0;i<n_city;i++){ //filling the first culumn
            X[i]=r*rnd.Rannyu(-1,1);
            Y[i]=r*rnd.Rannyu(-1,1);
            check+=i;
        }
    }
    if(in_circle==3){
        ifstream input("American_capitals.dat");
        string appo;
        if (input.is_open()){
            if(n_city!=50) cout<<"Not all American capitals"<<endl;
            input>>appo>>appo>>appo>>appo;
            for(int i=0;i<n_city;i++){
                input >>states[i]>>capitals[i]>> X[i]>>Y[i];
                check+=i;
            }
            input.close();
         } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    if(rank==0){
        cout <<"Position of cities"<<endl<<endl;
        cout <<"X"<<"           "<<"Y"<<endl;
        for(int i = 0; i< n_city; i++)
            cout <<X[i] << "    "<<Y[i]<<endl;
    }



    return ;}
void First_Travel(){
        for(int i=0;i<n_city;i++){ //filling matrix with int numbers
            for (int j=0;j<n_config;j++){
                A(i,j)=i;
            }
        }
   
    return;
}

arma::Mat<int> Mutation_1(arma::Mat<int> A){
    for (int j=1;j<n_config;j++){//swapping an element with the next in a culumn with 50% probability
        for(int i=2;i<n_city;i++){
           
            if(rnd.Rannyu()<0.3){
                swap(A(i,j),A(i-1,j));
            }
        }
    }
    return A;
}

arma::Mat<int> Mutation_2(arma::Mat<int> A){ //changing order of elements
    int m;
    for (int j=1;j<n_config;j++){
        m=int(rnd.Rannyu(3,n_city-1));
        int* values = new int[m];
        for(int i=n_city-m;i<n_city;i++){
            values[i-(n_city-m)]=A(i,j);
        }
        for(int i=n_city-m;i<n_city;i++){
            A(i,j)=values[n_city-i-1];
        }
        delete [] values;
    }
    return A;
}
   
void Check(){
    double a=A(0,0);
    for (int j=0;j<n_config;j++){
        if(A(0,j)!=a){
            cerr<<"Error: you start from the wrong city"<<endl;
            
        }
    }
    for(int j=0;j<n_config;j++){
        //double sum=accu(A.col(j));
        if(accu(A.col(j))!=check){ cerr<<"Error in swapping: you pass 2 times in the same city"<<endl;}
    }
    return ;
}



int Cost(){ //computing the distances for the culumn and select the min
    
    double* distance = new double[n_config];
    
    for (int j=0;j<n_config;j++){
        distance[j]=0;
        for(int i=0;i<n_city;i++){
            if(i<n_city-1) distance[j]+=Norm(X[A(i,j)],Y[A(i,j)],X[A(i+1,j)],Y[A(i+1,j)]);
            else distance[j]+=Norm(X[A(i,j)],Y[A(i,j)],X[A(0,j)],Y[A(0,j)]);
            
            }
    }
    
return Min(distance, n_config);
    delete[] distance;
}

double Norm(double Ax, double Ay, double Bx, double By){
    
    return fabs(Ay-By)+fabs(Ax-Bx);
}


template <class T> void swap ( T& a, T& b )
{
  T c(a); a=b; b=c;
}

int Min(double* a, int DIM){ //looking for min of an array and return index of culumn with shortest travel
    double min=a[0];
    int index=0;;
    for(int j=0;j<DIM;j++)
        {
            if(a[j]<min){
                min=a[j];
                index=j;
            }
        }
    return index;
}

void Compute_distance(int step,int j, bool print){ //calculatin the distance for a culumn
    dist=0;
    for(int i=0;i<n_city;i++){
        if(i<n_city-1) dist+=Norm(X[A(i,j)],Y[A(i,j)],X[A(i+1,j)],Y[A(i+1,j)]);
        else dist+=Norm(X[A(i,j)],Y[A(i,j)],X[A(0,j)],Y[A(0,j)]);
        }
    //if(dist<=best_dist) best_dist=dist;
    if(print){
        Distance<<step<<" "<<dist<<endl;
    }
    
    return;
}

void Reorganize_vec(int *vec){//copying vec in matrix every 2 culumn
    for (int j=0;j<n_config;j=j+2){
        for(int i=0;i<n_city;i++){
                A(i,j)=vec[i];
            }
    }
    
}
    
void Reorganize(){//copying the index_min culumn in matrix every 2 culumn
    for (int j=0;j<n_config;j=j+2){
        for(int i=0;i<n_city;i++){
            A(i,j)=A(i,index_min);
        }
    }
        
}


void Selection_operator(){
    i_ran=rnd.Rannyu(1,n_city-1);
    //j_ran=rnd.Rannyu(1,n_config-1); //don't change the first one
}

arma::Mat<int> Crossover(int i_ran,int j_ran){
    int l=0;
    int k=0;
    
    int* rna_j1_old = new int[n_city-i_ran];
    int* rna_j2_old = new int[n_city-i_ran];
    for(int i=i_ran;i<n_city;i++){
        rna_j1_old[i-i_ran]=A(i,j_ran);
        rna_j2_old[i-i_ran]=A(i,j_ran+1);
    }
    int* rna_j1_new = new int[n_city-i_ran];
    int* rna_j2_new = new int[n_city-i_ran];
    l=0;
    k=0;
    
    for(int i=1;i<n_city;i++){
        //cout<<search_in(A(5,j_ran+1),rna_j1_old)<<endl;
        if(search_in(A(i,j_ran+1),rna_j1_old)){ //if element of column is present in array save it in a new array in order
            rna_j1_new[l]=A(i,j_ran+1);
            l++;
        }
        if(search_in(A(i,j_ran),rna_j2_old)){ //the same in next column
            rna_j2_new[k]=A(i,j_ran);
            k++;
        }
    }
    for(int i=i_ran;i<n_city;i++){
        A(i,j_ran)=rna_j1_new[i-i_ran];
        A(i,j_ran+1)=rna_j2_new[i-i_ran];
    }
    delete [] rna_j1_old;
    delete [] rna_j1_new;
    delete [] rna_j2_old;
    delete [] rna_j2_new;
    return A;
}

bool search_in(int val,int *vec){
    bool a=false;
    for(int j=0; j<n_city-i_ran; j++){
        if (vec[j]-val==0)
            a=true; //it's present
    }
    
return a;
}

void Print_coordinates(){
    ofstream coordinates;
    if(in_circle==1) coordinates.open("output.cord_in_circle.dat",ios::app);
    if(in_circle==0) coordinates.open("output.coord_square.dat",ios::app);
    if(in_circle==3) coordinates.open("output.coord.dat",ios::app);
    for(int i=0;i<n_city;i++){
        coordinates<<states[A(i,0)]<<" "<<capitals[A(i,0)]<<" "<<X[A(i,0)]<<" "<<Y[A(i,0)]<<endl;
    }
    coordinates<<states[A(0,0)]<<" "<<capitals[A(0,0)]<<" "<<X[A(0,0)]<<" "<<Y[A(0,0)]<<endl;
    coordinates.close();
}

void Print_coordinates_city(int rank){
    for(int i=0;i<n_city;i++){
        coordinates<<X[A(i,0)]<<" "<<Y[A(i,0)]<<endl;
    }
    coordinates<<X[A(0,0)]<<" "<<Y[A(0,0)]<<endl;
}




