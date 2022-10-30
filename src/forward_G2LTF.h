// Forward modeling of vertical global-to-local (G2L) tansfer functions
// Consistent for magnetospheric and ionospheric Sq current systems

// References:
// Definition of G2L TFs: Tnm = -Bnm_{r}, see eq.(15) of:
// Püthe, C., Kuvshinov, A., & Olsen, N. (2015). Handling complex 
// source structures in global EM induction studies: From C-responses 
// to new arrays of transfer functions. Geophysical Journal 
// International, 201(1), 318–328. https://doi.org/10.1093/gji/ggv021

// Recursive formulation of Bnm_{r}, see Appendix H of :
// Alexey Kuvshinov, Alexey Semenov, Global 3-D imaging of mantle 
// electrical conductivity based on inversion of observatory 
// C-responses—I. An approach and its verification, Geophysical 
// Journal International, Volume 189, Issue 3, June 2012, Pages 
// 1335–1352, https://doi.org/10.1111/j.1365-246X.2011.05349.x

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _forward_G2LTF_H
#define _forward_G2LTF_H

#include "utils.h"
#include "legendre.h"
#include <vector>
#include <cmath>

inline void forward_G2LTF(std::vector<double> radius,
                          std::vector<double> conductivity,
                          std::vector<double> SH_degree,
                          std::vector<double> SH_order,
                          std::vector<double> periods,
                          double longitude, double colatitude,
                          std::vector<Dcomplex>& G2LTFs)
{
   // Loop all perioids
   G2LTFs.resize(periods.size());
   double theta = colatitude;
   double phi = longitude;
   double e_nm = 1.0; // 1nT unit source for G2L TFs
   for (int source_id=0; source_id<periods.size(); source_id++)
   {
      double n = SH_degree[source_id];
      double m = SH_order[source_id];
      double period = periods[source_id];
      double omega = 2*pi/period;

      /* Compute admittance Y1 of formula (E11) */
      // Compute intermediate parameters
      int N = radius.size();
      std::vector<Dcomplex> q(N);
      std::vector<Dcomplex> b(N);
      std::vector<Dcomplex> b_minus(N);
      std::vector<Dcomplex> b_plus(N);
      for (int k=0; k<N; k++)
      {
         q[k] = II*omega*mu0*radius[k];
         b[k] = std::sqrt((n+0.5)*(n+0.5)-II*omega*mu0
               *conductivity[k]*radius[k]*radius[k]);
         b_minus[k] = b[k]-0.5;
         b_plus[k] = b[k]+0.5;
      }
      std::vector<double> eta(N-1);
      std::vector<Dcomplex> zeta(N-1);
      std::vector<Dcomplex> tau(N-1);
      for (int k=0; k<N-1; k++)
      {
         eta[k] = radius[k]/radius[k+1];
         zeta[k] = std::pow(eta[k],2.0*b[k]);
         if (std::isinf(zeta[k].real())) tau[k] = Dcomplex(-1.0,0.);
         else tau[k] = (1.0-zeta[k])/(1.0+zeta[k]);
      }

      // Compute admittances Yl, formula (E11)
      std::vector<Dcomplex> Yl(N);
      Yl[N-1] = -b_plus[N-1]/q[N-1];
      for (int k=N-2; k>=0; k--)
      {
         Yl[k] = 1.0/q[k]*(q[k+1]*Yl[k+1]*(b[k]-0.5*tau[k]) 
               + b_plus[k]*b_minus[k]*tau[k])/((b[k]+0.5*tau[k]) 
               + q[k+1]*tau[k]*Yl[k+1]);
      }
      Dcomplex Y1 = Yl[0];

      /* Compute G2L TFs: Tnm = -Bnm_{r} */
      double x = cosd(theta);
      double Pnm = get_Pnm(n,std::abs(m),x);
      Dcomplex Ynm = Pnm*std::exp(II*m*(phi/180.0*pi));
      Dcomplex fac1 = II*omega*mu0*R0;
      Dcomplex fac2 = fac1*Y1;
      Dcomplex fac3 = fac2-n;
      Dcomplex Tnm = std::conj(-(2*n+1)*n*e_nm/fac3)*Ynm; // exp(-iwt) to exp(iwt)
      G2LTFs[source_id] = Tnm;
   }
}

#endif // _forward_G2LTF_H