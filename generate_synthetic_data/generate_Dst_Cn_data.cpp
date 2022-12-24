#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <random>
#include <iomanip>
#include <iomanip>
#include "../src/utils.h"
#include "../src/forward_Cn.h"

int main(int argc, char* argv[])
{
   if (argc < 2) 
   {
      std::cout << "Usage: " << argv[0] << " config_filename\n";
      return 1;
   }
   std::ifstream in_stream(argv[1]);
   assert(in_stream.good());
   std::string line;
   std::vector<std::string> para_str_vec;
   while (std::getline(in_stream, line))
   {
      if (*(line.begin())=='#' || line=="") continue;
      para_str_vec.push_back(line);
   }
   in_stream.close();

   double n;
   double period_min, period_max;
   int n_periods;
   double random_err;
   double data_uncertainty;
   int n_layers;
   std::stringstream ss;
   ss << para_str_vec[0]; ss >> n; ss.clear();
   ss << para_str_vec[1]; ss >> period_min; ss.clear();
   ss << para_str_vec[2]; ss >> period_max; ss.clear();
   ss << para_str_vec[3]; ss >> n_periods; ss.clear();
   ss << para_str_vec[4]; ss >> random_err; ss.clear();
   ss << para_str_vec[5]; ss >> data_uncertainty; ss.clear();
   ss << para_str_vec[6]; ss >> n_layers; ss.clear();
   std::vector<double> depth(n_layers);
   std::vector<double> radius(n_layers);
   std::vector<double> conductivity(n_layers);
   int idx = 7;
   for (int i=0; i<n_layers; i++)
   {
      ss << para_str_vec[idx]; 
      ss >> depth[i] >> conductivity[i]; 
      radius[i] = R0-depth[i];
      ss.clear();
      idx++;
   }
   std::ofstream synthetic_model_stream("synthetic_model.txt");
   for (int i=0; i<n_layers; i++)
   {
      synthetic_model_stream << depth[i] << "\t" << conductivity[i] << "\n";
   }
   synthetic_model_stream.close();

   // Generate log equally spaced periods
   std::vector<double> SH_degree(n_periods);
   std::vector<double> periods(n_periods);
   double a = log10(period_min);
   double b = log10(period_max);
   double step = (b-a)/(n_periods-1);
   for (int i=0; i<n_periods; i++)
   {
      SH_degree[i] = n;
      double period = a+i*step;
		periods[i] = pow(10, period);
   }

   // Forward modeling
   std::vector<Dcomplex> C_pred(n_periods);
   forward_Cn(radius, conductivity, SH_degree, periods, C_pred);

   // Add gaussian random noise
   std::vector<Dcomplex> C_pred_noise(n_periods);
   std::vector<double> d_err(n_periods);
   std::default_random_engine random_engine(time(0));
   std::normal_distribution<double> gaussian_distribution(0,1); // mu=0, std=1
   for (int i=0; i<n_periods; i++)
   {
      double rand = gaussian_distribution(random_engine);
      C_pred_noise[i] = C_pred[i]*(1+random_err*rand);
      d_err[i] = abs(C_pred[i])*data_uncertainty;
   }
   std::ofstream out_stream("synthetic_Dst_Cn_data.txt");
   for (int i=0; i<n_periods; i++)
   {
      out_stream << std::setw(4) << std::setiosflags(std::ios::fixed) << std::setprecision(0) << SH_degree[i]
                 << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << periods[i]
                 << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(6) << C_pred_noise[i].real()
                 << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(6) << C_pred_noise[i].imag()
                 << std::setw(15) << std::setiosflags(std::ios::fixed) << std::setprecision(6) << d_err[i] << "\n";
   }
   out_stream.close();

   return 0;
}