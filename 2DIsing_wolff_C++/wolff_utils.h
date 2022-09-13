#ifndef __WOLFF_UTILS__
#define __WOLFF_UTILS__
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <array>
#include<iomanip>

class wolff_utils
{
public: 
void left(int x, int y, int* left_x, int* left_y, int L)
{
    if (x < 1) 
    {
        *left_x = L-1; 
        *left_y = y; }
    else
    {
        *left_x = x-1;
        *left_y = y; 
    } 
} 

void right(int x, int y, int* right_x, int* right_y, int L)
{
    if (x > (L - 1.5)) 
    {
        *right_x = 0; 
        *right_y = y; }
    else
    {
        *right_x = x+1;
        *right_y = y; 
    } 

}

void up(int x, int y, int* up_x, int* up_y, int L)
{
    if (y < 0.5) 
    {
        *up_x = x; 
        *up_y = L-1; }
    else
    {
        *up_x = x;
        *up_y = y-1; 
    } 

}

void down(int x, int y, int* down_x, int* down_y, int L)
{
    if (y > (L-1.5)) 
    {
        *down_x = x; 
        *down_y = 0; }
    else
    {
        *down_x = x;
        *down_y = y+1; 
    } 
} 

}; 

#endif
