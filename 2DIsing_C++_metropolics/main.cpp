#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <array>
#include <iomanip>
#include <iostream>
#include "utils.h"
using namespace std;

utils UTILS; 

const int L = 40, B = 1, sweeps=50000, relax = 5000; //L is lattice size and B is Boltzmann constant;  
int M = 0, S = L*L; 
using Matrix = std::array<std::array<int, L>, L>;

//Inialize the L*L spin lattice 
Matrix InialSpin() {
    Matrix m = {};
    for(int i=0; i<L; i++)
    {
        for(int j=0; j<L; j++)
        {
            m[i][j] = 1; 
        }
    }
    return m;
}

int main()
{
    ofstream mfile;
    mfile.open("E:\\MonteCarlo\\2D-Ising-Model\\C++\\m1.dat");
    for (float T=0.5; T<4.0; T+=0.1)
    {
        //Inilazing the spin lattice SpinM 
        double Mag[sweeps+relax]; 
        double Mag_2[sweeps+relax]; 
        Matrix SpinM = InialSpin(); 
        for(int k=0; k<sweeps+relax; k++)
        {
            for(int i=0; i<(L*L); i++)
            {
                int tmp[8]; 
                for (int ww=0; ww<8; ww++)
                {
                    tmp[ww] = UTILS.randint(L);
                    //cout << UTILS.randint(L) << endl; 
                }
                int x_site = tmp[1];
                int y_site = tmp[4];
                //cout << x_site << endl;
                int top = SpinM[UTILS.Neigh1(x_site)][y_site];
                int bottom = SpinM[UTILS.Neigh2(x_site)][y_site];
                int left = SpinM[x_site][UTILS.Neigh1(y_site)];
                int right = SpinM[x_site][UTILS.Neigh2(y_site)];
                int delta_E = UTILS.funcdeltaE(SpinM[x_site][y_site], top, bottom, left, right); 
                if (delta_E <= 0) SpinM[x_site][y_site] = -SpinM[x_site][y_site]; 
                else if (delta_E > 0)
                {
                    if (exp(-delta_E/T) > UTILS.random()) SpinM[x_site][y_site] = -SpinM[x_site][y_site];
                }
            }
            int M = 0; 
            for(int i=0; i<L; i++)
            { for(int j=0; j<L; j++) M+=SpinM[i][j];}
            double m = (double)abs(M)/S;
            //double m2 = (double)abs(M)/S;
            Mag[k] = m; 
            Mag_2[k] = m*m; 
        }
        //Calculate the Magnetization 
        double M_avg = 0.0; //average of the spin lattice  
        double M_avg2 = 0.0; //
        double M_c = 0.0; 
        for(int w=relax; w<sweeps+relax; w++)
        {   
            M_avg+= Mag[w];
            M_avg2 += Mag_2[w]; 
        }
        M_avg = M_avg/sweeps;
        M_avg2 = M_avg2/sweeps;  
        M_c = (M_avg2 - M_avg*M_avg)/T;     
        mfile<<"   "<< T << "   "<< M_avg <<"   "<< M_c << endl; 
        cout << T << " " <<  M_avg << " " << M_c << endl; 
    }
    mfile.close();
    return 0; 
}