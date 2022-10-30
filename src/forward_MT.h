// Forward modeling of MT apparent resistivity and phase
// Modified from the matlab code of Digital Earth lab
// https://www.digitalearthlab.com/tutorial/tutorial-1d-mt-forward/

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _forward_MT_H
#define _forward_MT_H

#include "utils.h"
#include <cmath>
#include <vector>

inline void forward_MT(std::vector<double> thicknesses,
                       std::vector<double> resistivities,
                       std::vector<double> periods,
                       std::vector<double>& apparent_resistivity,
                       std::vector<double>& phase)
{
   apparent_resistivity.resize(periods.size());
   phase.resize(periods.size());
   int n = thicknesses.size();
   for (int f=0; f<periods.size(); f++)
   {
      double w = 2*pi/periods[f];
      std::vector<Dcomplex> impedances(n);
      Dcomplex Zn = std::sqrt(II*w*mu0*resistivities[n-1]);
      impedances[n-1] = Zn;
      for (int j=n-2; j>=0; j--)
      {
         double resistivity = resistivities[j];
         double thickness = thicknesses[j];
         Dcomplex dj = std::sqrt(II*(w*mu0*(1.0/resistivity)));
         Dcomplex wj = dj*resistivity;
         Dcomplex ej = std::exp(-2.0*thickness*dj);
         Dcomplex belowImpedance = impedances[j+1];
         Dcomplex rj = (wj-belowImpedance)/(wj+belowImpedance); 
         Dcomplex re = rj*ej;
         Dcomplex Zj = wj*((1.0-re)/(1.0+re));
         impedances[j] = Zj;
      }
      Dcomplex Z = impedances[0];
      double absZ = std::abs(Z);
      apparent_resistivity[f] = std::log10((absZ*absZ)/(w*mu0));
      phase[f] = std::atan2(std::imag(Z),std::real(Z))/pi*180.0;
   }
}

#endif // _forward_MT_H