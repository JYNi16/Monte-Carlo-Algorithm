#ifndef __UTILS__
#define __UTILS__
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <array>
#include<iomanip>
#define N 999
const int J = 1, l = 30; 
const float H = 0.0; 

class utils
{
public:
       int funcdeltaE(int loc, int top, int bottom, int left, int right)
       {
           int delta_E = 0; 
           int E_loc_former = -1 * loc * (J * (top + bottom + left + right) + H); 
           int E_loc_next =  loc * (J * (top + bottom + left + right) + H); 
           delta_E = E_loc_next - E_loc_former; 
           return delta_E; 
           }
       
       float random()
       {
           return rand()%(N+1)/(float)(N+1);
       }
       
       int randint(int Lattice)
       {
           return (rand()%(Lattice));
       }
       int Neigh1(int& Minus)
       {
           if ((Minus-1)%l == l-1 || (Minus-1)%l == -1) return (Minus-1+l);
           else return (Minus-1); 
       }

       int Neigh2(int& adds)
       {
           if( (adds+1)%l == 0 )  return (adds+1-l);
           else return (adds+1);
       }
};

#endif 