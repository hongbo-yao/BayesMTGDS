// Trans-dimensional Bayesian inversion of magnetotelluric and 
// geomagnetic depth sounding data

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#include <iostream>
#include "rjmcmc.h"
#include <mpi.h>

int main(int argc, char* argv[])
{
   MPI_Init(&argc, &argv);
   double tic = MPI_Wtime();
   int size, my_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if (argc < 2) 
   {
      if (my_rank==0) 
      {
         std::cout << "Usage: " << argv[0] << " config_filename\n";
      }
      MPI_Finalize();
      return 1;    
   }
   
   RJMCMC rjmcmc(argv[1]);
   rjmcmc.read_observed_data();
   rjmcmc.run_rjmcmc();
   MPI_Finalize();
   return 0;
}