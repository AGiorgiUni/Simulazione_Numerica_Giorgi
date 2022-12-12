#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"

#include "generatore.h"

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
    
    vector<double> sum(10000);
   
    int N=1;
    double appo=0;
    
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+rnd.Rannyu();
        }
        sum[i]=appo/N;
    }
    ofstream fout("Unif1.txt");
    
    for(int i=0;i<sum.size();i++){
        fout<<N<<" "<<sum[i]<<endl;
    }
    fout.close();
    
    ofstream fout_2("Unif2.txt");
    N=2;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+rnd.Rannyu();
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout_2<<N<<" "<<sum[i]<<endl;
    }
    fout_2.close();
    
    ofstream fout_10("Unif10.txt");
    N=10;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+rnd.Rannyu();
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout_10<<N<<" "<<sum[i]<<endl;
    }
    fout_10.close();
    
    ofstream fout_100("Unif100.txt");
    N=100;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+rnd.Rannyu();
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout_100<<N<<" "<<sum[i]<<endl;
    }
    fout_100.close();
    
    //Exponential distribution
    double lambda=1;
    N=1;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+ExpNumber(lambda,rnd.Rannyu());
        }
        sum[i]=appo/N;
    }
    ofstream fout2("Exp1.txt");
    
    for(int i=0;i<sum.size();i++){
        fout2<<N<<" "<<sum[i]<<endl;
    }
    fout2.close();
    
    ofstream fout2_2("Exp2.txt");
    
    N=2;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+ExpNumber(lambda,rnd.Rannyu());
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout2_2<<N<<" "<<sum[i]<<endl;
    }
    fout2_2.close();
    
    ofstream fout3_10("Exp10.txt");
    N=10;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+ExpNumber(lambda,rnd.Rannyu());
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout3_10<<N<<" "<<sum[i]<<endl;
    }
    fout3_10.close();
    
    ofstream fout4_100("Exp100.txt");
    N=100;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+ExpNumber(lambda,rnd.Rannyu());
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout4_100<<N<<" "<<sum[i]<<endl;
    }
    
    fout4_100.close();

    
    //Lorentzian Distribution
    double mu=0;
    double gamma=1;
    N=1;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+LorentzNumber(mu,rnd.Rannyu(),gamma);
        }
        sum[i]=appo/N;
    }
    ofstream fout3("Lorentz1.txt");
    
    for(int i=0;i<sum.size();i++){
        fout3<<N<<" "<<sum[i]<<endl;
    }
    fout3.close();
    
    ofstream fout5_2("Lorentz2.txt");
    N=2;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+LorentzNumber(mu,rnd.Rannyu(),gamma);
        }
        sum[i]=appo/N;
    }
    
    
    for(int i=0;i<sum.size();i++){
        fout5_2<<N<<" "<<sum[i]<<endl;
    }
    fout5_2.close();
    
    ofstream fout6_10("Lorentz10.txt");
    N=10;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+LorentzNumber(mu,rnd.Rannyu(),gamma);
        }
        sum[i]=appo/N;
    }
   
    
    for(int i=0;i<sum.size();i++){
        fout6_10<<N<<" "<<sum[i]<<endl;
    }
    
    fout6_10.close();
    
    ofstream fout7_100("Lorentz100.txt");
    N=100;
    for(int i=0;i<sum.size();i++){
        appo=0;
        for (int j=0;j<N;j++){
            appo=appo+LorentzNumber(mu,rnd.Rannyu(),gamma);
        }
        sum[i]=appo/N;
    }
  
    
    for(int i=0;i<sum.size();i++){
        fout7_100<<N<<" "<<sum[i]<<endl;
    }
    
    fout7_100.close();



    return 0;
}
