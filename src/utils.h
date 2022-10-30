// A collection of common variables and functions

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <complex>
#include <vector>
#include <stdexcept>
#include <cassert>

namespace EM 
{
   typedef std::complex<double>    Dcomplex;
   static const Dcomplex   II    = Dcomplex(0.,1.);
   static const double     pi    = 3.1415926535897932384626433832795;
   static const double     mu0   = 4.0*pi*1e-7;
   static const double     R0    = 6371200;  // Earth's mean radius in m
   static const double     z_cmb = 2890000;  // depth of core-mantle boundary in m
   static const double     sigma_core = 5e5; // conductivity of Earth's core, S/m
}
using namespace EM;

inline double cosd(double deg) 
{ 
   return cos(deg/180.0*pi); 
}


inline double sind(double deg) 
{
   return sin(deg/180.0*pi); 
}


template<class T>
inline int find_position(std::vector<T> vec, T number)
{
   // in default, vec is arranged from small to big
   int n_numbers = vec.size();
   if (number<vec[0] || number>vec[n_numbers-1])
   {
      // std::cout << "number: " << number << "\nvec:\n";
      // for (int i=0; i<n_numbers; i++)
      //    std::cout << vec[i] << "\n";
      throw std::invalid_argument("'number' must within 'vec'!");
   }
   int position = -1;
   for (int i=0; i<n_numbers-1; i++)
   {
      // if (number>=vec[i] && number<=vec[i+1])
      if (number>vec[i] && number<=vec[i+1])
      {
         position = i;
         break;
      }
   }
   assert(position!=-1);
   return position;
}


// Spherical harmonic sources for one period
struct SH_source
{
   double period;
   std::vector<double> SH_degree_list; // n
   std::vector<double> SH_order_list;  // m
   std::vector<Dcomplex> e_nm_list;
};


#endif // _UTILS_H