#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <array>
#include <iomanip>
#include <iostream>
#include <stack>
#include "utils.h"
#include "wolff_utils.h"
using namespace std;

utils UTILS; 
wolff_utils WOLFF_UTILS; 

const int L = 40, B = 1, sweeps=500, relax = 100; //L is lattice size and B is Boltzmann constant;  
int M = 0, S = L*L; 
using Matrix = std::array<std::array<int, L>, L>;
bool IsOnes = true; 
stack <int> ls; 

//Inialize the L*L spin lattice 
Matrix InialSpin(bool IsOnes) {
    Matrix m = {};
    for(int i=0; i<L; i++)
    {
        for(int j=0; j<L; j++)
        {
            if (IsOnes) {m[i][j] = 1;}
            else 
            {
                if (UTILS.random() < 0.5)
                {m[i][j] = 1;}
                else {m[i][j] = -1;}
            } 
        }
    }
    return m;
}

int main()
{
    ofstream mfile;
    mfile.open("E:\\MonteCarlo\\2D-Ising-Model\\C++\\m1.dat");
    int left_x, left_y; 
    int right_x, right_y; 
    int up_x, up_y; 
    int down_x, down_y; 

    for (float T=0.5; T<4.0; T+=0.1)
    {
        //Inilazing the spin lattice SpinM 
        float P_add = (1 - exp(-2 * J / T));
        double Mag[sweeps]; 
        double Mag_2[sweeps]; 
        Matrix SpinM = InialSpin(!IsOnes); 
        for(int k=0; k<sweeps+relax; k++)
        {
            int x = UTILS.randint(L);
            int y = UTILS.randint(L);
            int sign = SpinM[x][y]; 
            ls.push(x);
            ls.push(y); 
            Matrix lable = InialSpin(! IsOnes);
            lable[x][y] = 0;

            while(! ls.empty())
            {
                int y_site = ls.top(); 
                //cout << "y_site is:" << y_site << endl; 
                ls.pop(); 
                int x_site = ls.top(); 
                ls.pop(); 
                SpinM[x_site][y_site] = -sign;
                
                //left
                WOLFF_UTILS.left(x_site, y_site, &left_x, &left_y, L); 
                if (SpinM[left_x][left_y] * sign > 0.5 && UTILS.random() < P_add)
                {
                    ls.push(left_x); 
                    ls.push(left_y); 
                    lable[left_x][left_y] =0; 
                }

                //right 
                WOLFF_UTILS.right(x_site, y_site, &right_x, &right_y, L); 
                if (SpinM[right_x][right_y] * sign > 0.5 && UTILS.random() < P_add)
                {
                    ls.push(right_x); 
                    ls.push(right_y); 
                    lable[right_x][right_y] =0; 
                }

                //up 
                WOLFF_UTILS.up(x_site, y_site, &up_x, &up_y, L); 
                if (SpinM[up_x][up_y] * sign > 0.5 && UTILS.random() < P_add)
                {
                    ls.push(up_x); 
                    ls.push(up_y); 
                    lable[up_x][up_y] =0; 
                }

                //down 
                WOLFF_UTILS.down(x_site, y_site, &down_x, &down_y, L); 
                if (SpinM[down_x][down_y] * sign > 0.5 && UTILS.random() < P_add)
                {
                    ls.push(down_x); 
                    ls.push(down_y); 
                    lable[down_x][down_y] =0; 
                }
            }

            if (k>=relax)
            {
                int M = 0;
                for(int i=0; i<L; i++)
                { 
                    for(int j=0; j<L; j++) 
                    M+=SpinM[i][j];
                }

                //cout << "M is" << M << endl; 
                double m = (double)abs(M)/S;
                Mag[k-relax] = m;
                Mag_2[k-relax] = m*m; 
            }
        }
        //cout << "Mag_2 is:" << Mag << endl; 
        //Calculate the Magnetization 
        double M_avg = 0.0; //average of the spin lattice  
        double M_avg2 = 0.0; //
        double M_c = 0.0; 
        for (int w=0; w<sweeps; w++)
        {
            M_avg += Mag[w]; 
            M_avg2 += Mag_2[w]; 
        }
        M_avg /= sweeps;
        M_avg2 /= sweeps;  
        M_c = (M_avg2 - M_avg*M_avg)/T;     
        mfile<<"   "<< T << "   "<< M_avg <<"   "<< M_c << endl; 
        cout << T << " " <<  M_avg << " " << M_c << endl; 
    }
    mfile.close();
    return 0; 
}