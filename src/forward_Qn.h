// Forward modeling of geomagnetic Qn-responses 
// for a spherically layered conductivity model

// Reference:
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

#ifndef _forward_Qn_H
#define _forward_Qn_H

#include "utils.h"
#include <cmath>
#include <vector>

inline void forward_Qn(std::vector<double> radius,
                       std::vector<double> conductivity,
                       std::vector<double> SH_degree,
                       std::vector<double> periods,
                       std::vector<Dcomplex>& Qn_responses)
{
   int M = radius.size();
   Qn_responses.resize(periods.size());
   for (int f=0; f<periods.size(); f++)
   {
      double n = SH_degree[f];
      double omega = 2*pi/periods[f];

      // Compute intermediate parameters
      std::vector<Dcomplex> q(M);
      std::vector<Dcomplex> b(M);
      std::vector<Dcomplex> b_minus(M);
      std::vector<Dcomplex> b_plus(M);
      for (int k=0; k<M; k++)
      {
         q[k] = II*omega*mu0*radius[k];
         b[k] = std::sqrt((n+0.5)*(n+0.5)-II*omega*mu0
               *conductivity[k]*radius[k]*radius[k]);
         b_minus[k] = b[k]-0.5;
         b_plus[k] = b[k]+0.5;
      }
      std::vector<double> eta(M-1);
      std::vector<Dcomplex> zeta(M-1);
      std::vector<Dcomplex> tau(M-1);
      for (int k=0; k<M-1; k++)
      {
         eta[k] = radius[k]/radius[k+1];
         zeta[k] = std::pow(eta[k],2.0*b[k]);
         if (std::isinf(zeta[k].real())) tau[k] = Dcomplex(-1.0,0.);
         else tau[k] = (1.0-zeta[k])/(1.0+zeta[k]);
      }

      // Compute admittances Yl, formula (E11)
      std::vector<Dcomplex> Yl(M);
      Yl[M-1] = -b_plus[M-1]/q[M-1];
      for (int k=M-2; k>=0; k--)
      {
         Yl[k] = 1.0/q[k]*(q[k+1]*Yl[k+1]*(b[k]-0.5*tau[k]) 
               + b_plus[k]*b_minus[k]*tau[k])/((b[k]+0.5*tau[k]) 
               + q[k+1]*tau[k]*Yl[k+1]);
      }
      Dcomplex Y1 = Yl[0];

      // Compute Qn-responses
      Dcomplex Qn = n/(n+1.0)*(II*omega*mu0*R0*Y1+n+1.0)/(II*omega*mu0*R0*Y1-n);
      Qn = std::conj(Qn); // convert from exp(-iwt) to exp(iwt)
      Qn_responses[f] = Qn;
   }
}

#endif // _forward_Qn_H
