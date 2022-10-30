// Compute complex spherical harmonics and gradient based on Schmidt 
// quasi-normalized associated Legendre function of degree n and order |m|

// Y_{n}^{m}(\theta,\phi) = P_{n}^{|m|}(\cos\theta) \exp(im\phi)
// Has been checked by comparing with Matlab

// Reference: e.g.
// Sabaka T.J., Hulot G., Olsen N. (2013) Mathematical Properties Relevant 
// to Geomagnetic Field Modeling. In: Freeden W., Nashed M., Sonar T. (eds) 
// Handbook of Geomathematics. Springer, Berlin, Heidelberg. 
// https://doi.org/10.1007/978-3-642-27793-1_17-2

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _LEGENDRE_H
#define _LEGENDRE_H

#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>

inline double Rnm(double n, double m)
{
   return sqrt(n*n-m*m);
}



inline double get_Pnm(double n, double m, double x)
{
   // x=cos(theta)
   // s=sin(theta)=sqrt(1-x^2);
   double s = sqrt(1-x*x);
   if (n==0 && m==0) return 1.0;
   else if (n==1 && m==0) return x;
   else if (n==1 && m==1) return s;
   else if (n==m && n>1) 
   {
      return sqrt((2*n-1)/(2*n))*s*get_Pnm(n-1,n-1,x);
   }
   else if (n>m && m>=0)
   {
      if (n-1==m)
         return (2*n-1)*x*get_Pnm(n-1,m,x)/Rnm(n,m);
      else
         return ((2*n-1)*x*get_Pnm(n-1,m,x)-Rnm(n-1,m)*get_Pnm(n-2,m,x))/Rnm(n,m);
   }
   else std::abort();
}



inline double get_dPnm(double n, double m, double x)
{
   double s = sqrt(1-x*x);
   return (n*x*get_Pnm(n,m,x)-Rnm(n,m)*get_Pnm(n-1,m,x))/s;
}



// Compute Pnm and dPnm simultaneously
inline std::vector<double> legendre(double n, double m, double x)
{
   double Pnm = 0.0;
   double dPnm = 0.0;

   // x=cos(theta)
   // s=sin(theta)=sqrt(1-x^2);
   double s = sqrt(1-x*x);
   if (n==0 && m==0)
   {
      Pnm = 1.0;
      dPnm = 0.0;
   }
   else if (n==1 && m==0)
   {
      Pnm = x;
      dPnm = -s;
   }
   else if (n==1 && m==1)
   {
      Pnm = s;
      dPnm = x;
   }
   else if (n==m && n>1)
   {
      Pnm = sqrt((2*n-1)/(2*n))*s*legendre(n-1,n-1,x)[0];
      dPnm = (n*x*Pnm)/s;
   }
   else if (n>m && m>=0)
   {
      if (n-1==m)
         Pnm = (2*n-1)*x*legendre(n-1,m,x)[0]/Rnm(n,m);
      else
         Pnm = ( (2*n-1)*x*legendre(n-1,m,x)[0]-Rnm(n-1,m)*legendre(n-2,m,x)[0] )/Rnm(n,m);
      dPnm = ( n*x*Pnm-Rnm(n,m)*legendre(n-1,m,x)[0] )/s;
   }
   std::vector<double> value(2);
   value[0] = Pnm;
   value[1] = dPnm;
   return value;
}

#endif // _LEGENDRE_H