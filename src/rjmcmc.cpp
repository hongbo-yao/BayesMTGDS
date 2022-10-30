// Trans-dimensional Bayesian inversion of magnetotelluric and 
// geomagnetic depth sounding data

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#include "rjmcmc.h"
#include "forward_Qn.h"
#include "forward_Cn.h"
#include "forward_G2LTF.h"
#include "forward_MT.h"

#include <sys/stat.h> // int mkdir(const char*, __mode_t)
#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <mpi.h>

// C++ 11 random number
#include <random>
#include <ctime>

std::random_device rd;
std::default_random_engine random_engine(rd());
std::uniform_real_distribution<double> uniform_distribution(0,1);
std::normal_distribution<double> gaussian_distribution(0,1); // mu=0, std=1

RJMCMC::RJMCMC(std::string config_file)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   read_input_parameters(config_file);
   if (my_rank==0) print_computation_parameters();
}



RJMCMC::~RJMCMC()
{
   
}



void RJMCMC::read_input_parameters(std::string config_file)
{
   std::ifstream in_stream(config_file);
   assert(in_stream.good());
   std::string line;
   std::string key,value,comment;
   bool no_given_h_min = true;
   bool no_given_birth_death_std = true;
   bool no_given_update_std = true;
   bool no_given_move_std = true;
   while (std::getline(in_stream, line))
   {
      // skip line that starts with '#' or empty
      if (*(line.begin()) == '#' || line == "") continue;
      int loc = line.find(':');
      if(loc != std::string::npos)
      {
         // obtain key
         key = line.substr(0,loc);
         std::transform(key.begin(),key.end(),key.begin(),tolower);

         // obtain value
         int comment_loc = line.find('#');
         if(comment_loc == std::string::npos) // no comments
            value = line.substr(loc+1,line.size());
         else
            value = line.substr(loc+1,comment_loc);

         // set value
         if (key.find("observed data file")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> observed_data_file;
         }
         if (key.find("minimum number of layers")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> k_min;
         }
         if (key.find("maximum number of layers")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> k_max;
         }
         if (key.find("minimum interface depth")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> z_min;
         }
         if (key.find("maximum interface depth")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> z_max;
         }
         if (key.find("minimum log10 conductivity")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> sigma_min;
         }
         if (key.find("maximum log10 conductivity")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> sigma_max;
         }
         if (key.find("bitrh probability")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> birth_prob;
         }
         if (key.find("death probability")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> death_prob;
         }
         if (key.find("move probability")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> move_prob;
         }
         if (key.find("update probability")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> update_prob;
         }
         if (key.find("number of mcmc iterations")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> maxit;
         }
         if (key.find("print every iterations")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> print_every_iteration;
         }
         if (key.find("minimun layer thickness")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> h_min;
            no_given_h_min = false;
         }
         if (key.find("prior distribution type")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> prior_type;
         }
         if (key.find("smooth factor")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> log10_conductivity_gradient_std;
         }
         if (key.find("birth/death pertubation")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> birth_death_std;
            no_given_birth_death_std = false;
         }
         if (key.find("update pertubation")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> update_std;
            no_given_update_std = false;
         }
         if (key.find("move pertubation")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> move_std;
            no_given_move_std = false;
         }
         if (key.find("burn in period")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> burnin;
         }
         if (key.find("thin samples")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> thin;
         }
         if (key.find("depth resolution")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_depth_blocks;
         }
         if (key.find("conductivity resolution")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_sigma_blocks;
         }
         if (key.find("credibl interval")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> credible_interval;
         }
         if (key.find("model interface depth")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> model_interface_depth_file;
         }
         if (key.find("coarse depth number")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_depth_blocks_coarse;
         }
         if (key.find("write predicted responses")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> write_predicted_responses;
         }
         if (key.find("write rms misfit")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> write_rms_misfit;
         }
         if (key.find("output directory")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> output_dir;
         }
         if (key.find("write model samples")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> write_model_samples;
         }
      }
   }

   // Compute default parameters
   if (no_given_h_min) h_min = (z_max-z_min)/(2*k_max);
   if (no_given_birth_death_std) birth_death_std = (sigma_max-sigma_min)*0.1;
   if (no_given_update_std) update_std = (sigma_max-sigma_min)*0.15;
   if (no_given_move_std) move_std = (z_max-z_min)*0.1;
   left_credible_interval = 1-credible_interval;
   if (prior_type==1)
   {
      bestfit_halfspace_model = get_bestfit_halfspace_model();
   }

   // Check parameters
   if (observed_data_file=="" && my_rank==0)
   {
      throw std::invalid_argument("Error: no observed data file!\n");
   }
   assert(k_min>=1 && k_min<k_max);
   assert(z_min<z_max && z_max<=mcmc_z_cmb);
   assert(sigma_min<sigma_max);
   delta_sigma = sigma_max-sigma_min;
   double prob_sum = birth_prob+death_prob+move_prob+update_prob;
   if (std::abs(prob_sum-1)>1e-12)
   {
      if (my_rank==0) std::cout << "Warning: prob_sum!=1.0, default values have been used\n";
      birth_prob = 0.25;
      death_prob = 0.25;
      move_prob = 0.25;
      update_prob = 0.25;
   }
   if (birth_prob!=death_prob)
   {
      if (my_rank==0) std::cout << "Warning: birth_prob!=death_prob, default values have been used\n";
      birth_prob = 0.25;
      death_prob = 0.25;
      move_prob = 0.25;
      update_prob = 0.25;
   }
   assert(h_min>0);
   assert(maxit>print_every_iteration);
   assert(log10_conductivity_reference_std>0);
   assert(log10_conductivity_gradient_std>0);
   assert(birth_death_std>0);
   assert(update_std>0);
   assert(move_std>0);
   assert(burnin<maxit);

   // Make output directory and print plotting parameters
   if (my_rank==0)
   {
      mkdir(output_dir.c_str(),S_IRWXU);
      std::string out_filename = output_dir + "/plot_parameters.dat";
      std::ofstream out_stream(out_filename);
      out_stream << z_min << "\n";
      out_stream << z_max << "\n";
      out_stream << maxit << "\n";
      out_stream << burnin << "\n";
      out_stream.close();
   }

   // Set save flag
   if (write_predicted_responses==1) // only write responses of rank 0
   {
      if (my_rank==0) save_predicted_responses = true;
   }
   else if (write_predicted_responses==2) // write responses of all processes
   {
      save_predicted_responses = true;
   }

   if (write_rms_misfit==1) // only write rms of rank 0
   {
      if (my_rank==0) save_rms_misfit = true;
   }
   else if (write_rms_misfit==2) // write rms of all processes
   {
      save_rms_misfit = true;
   }

   // check file stream
   if (model_interface_depth_file!="" && my_rank==0)
   {
      std::ifstream model_in_stream(model_interface_depth_file);
      if (!model_in_stream.good())
         throw std::invalid_argument("Invalid model_interface_depth_file\n");
      model_in_stream.close();
   }
}



void RJMCMC::print_computation_parameters(std::string out)
{
   std::ofstream out_stream(out);
   out_stream <<"# -------------------------\n";
   out_stream <<"# Observed data:\n";
   out_stream <<"# -------------------------\n";
   out_stream << "Observed data file            :  " << observed_data_file << "\n";
   
   out_stream << "\n";
   out_stream <<"# -------------------------\n";
   out_stream <<"# Common RJMCMC parameters:\n";
   out_stream <<"# -------------------------\n";
   out_stream << "Minimum number of layers      :  " << k_min << "\n";
   out_stream << "Maximum number of layers      :  " << k_max << "\n";
   out_stream << "Minimum interface depth       :  " << z_min << "\n";
   out_stream << "Maximum interface depth       :  " << z_max << "\n";
   out_stream << "Minimum log10 conductivity    :  " << sigma_min << "\n";
   out_stream << "Maximum log10 conductivity    :  " << sigma_max << "\n";
   out_stream << "Bitrh probability             :  " << birth_prob << "\n";
   out_stream << "Death probability             :  " << death_prob << "\n";
   out_stream << "Move probability              :  " << move_prob << "\n";
   out_stream << "Update probability            :  " << update_prob << "\n";
   out_stream << "Number of mcmc iterations     :  " << maxit << "\n";
   out_stream << "Print every iterations        :  " << print_every_iteration << "\n";

   out_stream << "\n";
   out_stream <<"# -------------------------\n";
   out_stream <<"# Prior parameters:\n";
   out_stream <<"# -------------------------\n";
   out_stream << "Minimun layer thickness       :  " << h_min << "\n";
   out_stream << "Prior distribution type       :  " << prior_type << "\n";
   out_stream << "Smooth factor                 :  " << log10_conductivity_gradient_std << "\n";

   out_stream << "\n";
   out_stream <<"# -------------------------\n";
   out_stream <<"# Proposal parameters:\n";
   out_stream <<"# -------------------------\n";
   out_stream << "Birth/death pertubation       :  " << birth_death_std << "\n";
   out_stream << "Update pertubation            :  " << update_std << "\n";
   out_stream << "Move pertubation              :  " << move_std << "\n";

   out_stream << "\n";
   out_stream <<"# -------------------------\n";
   out_stream <<"# Post-processing:\n";
   out_stream <<"# -------------------------\n";
   out_stream << "Burn in period                :  " << burnin << "\n";
   out_stream << "Thin samples                  :  " << thin << "\n";
   out_stream << "Depth resolution              :  " << n_depth_blocks << "\n";
   out_stream << "Conductivity resolution       :  " << n_sigma_blocks << "\n";
   out_stream << "Credibl interval              :  " << credible_interval << "\n";
   if (model_interface_depth_file!="")
      out_stream << "Model interface depth         :  " << model_interface_depth_file << "\n";
   else
      out_stream << "Coarse depth resolution       :  " << n_depth_blocks_coarse << "\n";
   out_stream << "Write predicted responses     :  " << write_predicted_responses << "\n";
   out_stream << "Write rms misfit              :  " << write_rms_misfit << "\n";
   out_stream << "Output directory              :  " << output_dir << "\n";
   out_stream << "Write model samples           :  " << write_model_samples << "\n";

   out_stream.close();
}



void RJMCMC::read_observed_data()
{
      // Clear old data
   station_name = "";
   GG_longitude = 9999;
   GG_colatitude = 9999;
   GM_longitude = 9999;
   GM_colatitude = 9999;

   N_data = 0;
   MT_period_info.clear();
   MT_all_periods.clear();
   MT_appres_Ndata = 0;
   MT_appres_period_id.clear();
   MT_appres_periods.clear();
   MT_appres.clear();
   MT_appres_uncertainty.clear();
   MT_phase_Ndata = 0;
   MT_phase_period_id.clear();
   MT_phase_periods.clear();
   MT_phase.clear();
   MT_phase_uncertainty.clear();

   Qn_Ndata = 0;
   Qn_periods.clear();
   Qn_SH_degree.clear();
   Qn_responses.clear();
   Qn_data_uncertainty.clear();

   Cn_Ndata = 0;
   Cn_periods.clear();
   Cn_SH_degree.clear();
   Cn_responses.clear();
   Cn_data_uncertainty.clear();

   Dst_G2LTF_Ndata = 0;
   Dst_G2LTF_periods.clear();
   Dst_G2LTF_SH_degree.clear();
   Dst_G2LTF_SH_order.clear();
   Dst_G2LTF_responses.clear();
   Dst_G2LTF_data_uncertainty.clear();

   Sq_G2LTF_Ndata = 0;
   Sq_G2LTF_periods.clear();
   Sq_G2LTF_SH_degree.clear();
   Sq_G2LTF_SH_order.clear();
   Sq_G2LTF_responses.clear();
   Sq_G2LTF_data_uncertainty.clear();

   // Read new data
   if (my_rank==0)
   {
      std::cout << "   ========================================== Data ==========================================\n";
      std::cout << "   # TF_type  period_id    period       n     m         real           imag         std_err\n";
   }
   std::ifstream in_stream(observed_data_file);
   assert(in_stream.good());
   std::string line;
   std::string key,value,comment;
   int n_data_readin = 0;
   while (std::getline(in_stream, line))
   {
      // skip line that starts with '#' or empty
      if (*(line.begin()) == '#' || line == "") continue;
      int loc = line.find(':');
      if(loc != std::string::npos)
      {
         // obtain key
         key = line.substr(0,loc);
         std::transform(key.begin(),key.end(),key.begin(),tolower);

         // obtain value
         int comment_loc = line.find('#');
         if(comment_loc == std::string::npos) // no comments
            value = line.substr(loc+1,line.size());
         else
            value = line.substr(loc+1,comment_loc);

         // set value
         if (key.find("station name")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> station_name;
         }
         if (key.find("gg longitude")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> GG_longitude;
         }
         if (key.find("gg latitude")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            double latitude;
            value_sstr >> latitude;
            GG_colatitude = 90-latitude;
         }
         if (key.find("gm longitude")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> GM_longitude;
         }
         if (key.find("gm latitude")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            double latitude;
            value_sstr >> latitude;
            GM_colatitude = 90-latitude;
         }
         if (key.find("number of data")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> N_data;
         }
      }

      if (loc==std::string::npos && N_data>0 && n_data_readin<N_data)
      {
         if (my_rank==0)
         {
            std::cout << "   " << line << "\n";
         }
         std::stringstream sstr(line);
         std::string TF_type;
         int period_id;
         double period, n, m, TF_re, TF_im, TF_std_err;
         sstr >> TF_type >> period_id >> period >> n >> m >> TF_re >> TF_im >> TF_std_err;
         if (TF_type=="Rho") // MT apparent resistivity
         {
            MT_period_info.insert(std::make_pair(period_id,period));
            MT_appres_period_id.push_back(period_id);
            MT_appres_periods.push_back(period);
            MT_appres.push_back(TF_re);
            MT_appres_uncertainty.push_back(TF_std_err);
         }
         else if (TF_type=="Phase") // MT impedance phase
         {
            MT_period_info.insert(std::make_pair(period_id,period));
            MT_phase_period_id.push_back(period_id);
            MT_phase_periods.push_back(period);
            MT_phase.push_back(TF_re);
            MT_phase_uncertainty.push_back(TF_std_err);
         }
         else if (TF_type=="Q") // global Q-responses
         {
            Qn_periods.push_back(period);
            Qn_SH_degree.push_back(n);
            Qn_responses.push_back(Dcomplex(TF_re,TF_im));
            Qn_data_uncertainty.push_back(TF_std_err);
         }
         else if (TF_type=="C") // local C-responses
         {
            Cn_periods.push_back(period);
            Cn_SH_degree.push_back(n);
            Cn_responses.push_back(Dcomplex(TF_re,TF_im));
            Cn_data_uncertainty.push_back(TF_std_err);
         }
         else if (TF_type=="SqG2L") // Sq global-to-local TFs
         {
            Sq_G2LTF_periods.push_back(period);
            Sq_G2LTF_SH_degree.push_back(n);
            Sq_G2LTF_SH_order.push_back(m);
            Sq_G2LTF_responses.push_back(Dcomplex(TF_re,TF_im));
            Sq_G2LTF_data_uncertainty.push_back(TF_std_err);
         }
         else if (TF_type=="DstG2L") // Dst global-to-local TFs
         {
            Dst_G2LTF_periods.push_back(period);
            Dst_G2LTF_SH_degree.push_back(n);
            Dst_G2LTF_SH_order.push_back(m);
            Dst_G2LTF_responses.push_back(Dcomplex(TF_re,TF_im));
            Dst_G2LTF_data_uncertainty.push_back(TF_std_err);
         }
         n_data_readin++;
      }
   }

   MT_appres_Ndata = MT_appres_period_id.size();
   MT_phase_Ndata = MT_phase_period_id.size();
   Qn_Ndata = Qn_periods.size();
   Cn_Ndata = Cn_periods.size();
   Dst_G2LTF_Ndata = Dst_G2LTF_periods.size();
   Sq_G2LTF_Ndata = Sq_G2LTF_periods.size();

   if (Dst_G2LTF_Ndata>0 && (GM_longitude==9999||GM_colatitude==9999) && my_rank==0)
   {
      throw std::invalid_argument("Inverting Dst G2LTF needs geomagnetic longitude and latitude!\n");
   }

   if (Sq_G2LTF_Ndata>0 && (GG_longitude==9999||GG_colatitude==9999) && my_rank==0)
   {
      throw std::invalid_argument("Inverting Sq G2LTF needs geographic longitude and latitude!\n");
   }

   // Get all MT periods
   for (std::map<int,double>::iterator i=MT_period_info.begin(); i!=MT_period_info.end(); i++)
   {
      MT_all_periods.push_back(i->second);
   }

   // Write observed data for plotting
   if (MT_appres_Ndata>0 && my_rank==0)
   {
      std::string out_filename = output_dir + "/observed_MT_appres_responses.dat";
      std::ofstream observed_out_stream(out_filename);
      for (int i=0; i<MT_appres_Ndata; i++)
      {
         observed_out_stream  << MT_appres_periods[i] << "\t" 
                              << MT_appres[i] << "\t" 
                              << MT_appres_uncertainty[i] << "\n";
      }
      observed_out_stream.close();
   }

   if (MT_phase_Ndata>0 && my_rank==0)
   {
      std::string out_filename = output_dir + "/observed_MT_phase_responses.dat";
      std::ofstream observed_out_stream(out_filename);
      for (int i=0; i<MT_phase_Ndata; i++)
      {
         observed_out_stream  << MT_phase_periods[i] << "\t" 
                              << MT_phase[i] << "\t" 
                              << MT_phase_uncertainty[i] << "\n";
      }
      observed_out_stream.close();
   }

   if (Qn_Ndata>0 && my_rank==0)
   {
      std::string out_filename = output_dir + "/observed_Qn_responses.dat";
      std::ofstream observed_out_stream(out_filename);
      for (int i=0; i<Qn_Ndata; i++)
      {
         observed_out_stream  << Qn_SH_degree[i] << "\t"
                              << Qn_periods[i] << "\t" 
                              << Qn_responses[i].real() << "\t" 
                              << Qn_responses[i].imag() << "\t" 
                              << Qn_data_uncertainty[i] << "\n";
      }
      observed_out_stream.close();
   }

   if (Cn_Ndata>0 && my_rank==0)
   {
      std::string out_filename = output_dir + "/observed_Cn_responses.dat";
      std::ofstream observed_out_stream(out_filename);
      for (int i=0; i<Cn_Ndata; i++)
      {
         observed_out_stream  << Cn_SH_degree[i] << "\t"
                              << Cn_periods[i] << "\t" 
                              << Cn_responses[i].real() << "\t" 
                              << Cn_responses[i].imag() << "\t" 
                              << Cn_data_uncertainty[i] << "\n";
      }
      observed_out_stream.close();
   }

   if (Dst_G2LTF_Ndata>0 && my_rank==0)
   {
      std::string out_filename = output_dir + "/observed_Dst_G2LTF_responses.dat";
      std::ofstream observed_out_stream(out_filename);
      for (int i=0; i<Dst_G2LTF_Ndata; i++)
      {
         observed_out_stream  << Dst_G2LTF_SH_degree[i] << "\t"
                              << Dst_G2LTF_SH_order[i] << "\t"
                              << Dst_G2LTF_periods[i] << "\t" 
                              << Dst_G2LTF_responses[i].real() << "\t" 
                              << Dst_G2LTF_responses[i].imag() << "\t" 
                              << Dst_G2LTF_data_uncertainty[i] << "\n";
      }
      observed_out_stream.close();
   }

   if (Sq_G2LTF_Ndata>0 && my_rank==0)
   {
      std::string out_filename = output_dir + "/observed_Sq_G2LTF_responses.dat";
      std::ofstream observed_out_stream(out_filename);
      for (int i=0; i<Sq_G2LTF_Ndata; i++)
      {
         observed_out_stream  << Sq_G2LTF_SH_degree[i] << "\t"
                              << Sq_G2LTF_SH_order[i] << "\t"
                              << Sq_G2LTF_periods[i] << "\t" 
                              << Sq_G2LTF_responses[i].real() << "\t" 
                              << Sq_G2LTF_responses[i].imag() << "\t" 
                              << Sq_G2LTF_data_uncertainty[i] << "\n";
      }
      observed_out_stream.close();
   }
}



double RJMCMC::get_bestfit_halfspace_model()
{
    // linspace of conductivity
   std::vector<double> sigma_linspace(n_sigma_blocks+1);
   double sigma_interval = delta_sigma/n_sigma_blocks;
   for (int i=0; i<sigma_linspace.size(); i++)
   {
      sigma_linspace[i] = sigma_min + i*sigma_interval;
   }

   std::vector<double> misfit(n_sigma_blocks+1);
   for (int i=0; i<misfit.size(); i++)
   {
      std::vector<double> z;
      std::vector<double> sigma(1,sigma_linspace[i]);
      struct PredictedResponse temp_predicted_response;
      double temp_misfit;
      double temp_rms;
      forward_modeling(z, sigma, temp_predicted_response, temp_misfit, temp_rms);
      misfit[i] = temp_misfit;
   }

   auto min_pdf_pos = std::min_element(misfit.begin(), misfit.end());
   return sigma_linspace[min_pdf_pos-misfit.begin()];
}



double RJMCMC::get_prior_probability(std::vector<double> z, std::vector<double> sigma)
{
   int k = sigma.size();

   if (prior_type==0) // uniform distribution
   {
      return 1.0/std::pow(delta_sigma,k);
   }

   else if (prior_type==1) // normal distribution with reference model constraint
   {
      double prior_covariance = log10_conductivity_reference_std*log10_conductivity_reference_std;
      double p1 = std::pow(2*pi,k)*std::pow(prior_covariance,k);
      p1 = std::pow(p1,-0.5);
      double sum = 0.;
      for (int i=0; i<k; i++)
      {
         sum += std::pow(sigma[i]-bestfit_halfspace_model,2.0)/prior_covariance;
      }
      double p2 = std::exp(-0.5*sum);
      return p1*p2;
   }

   else // normal distribution with smooth constraint
   {
      double prior_covariance = log10_conductivity_gradient_std*log10_conductivity_gradient_std;
      double p1 = std::pow(2*pi,k-1)*std::pow(prior_covariance,k-1);
      p1 = std::pow(p1,-0.5);

      // std::vector<double> h(k-1);
      // h[0] = z[0];
      // for (int i=1; i<k-1; i++)
      // {
      //    h[i] = z[i]-z[i-1];
      // }
      // double sum = 0.;
      // for (int i=0; i<k-1; i++)
      // {
      //    double dsigma_dz = (sigma[i+1]-sigma[i])/(h[i]-h_min);
      //    sum += dsigma_dz*dsigma_dz/prior_covariance;
      // }

      double sum = 0.;
      for (int i=0; i<k-1; i++)
      {
         sum += std::pow(sigma[i+1]-sigma[i],2.0)/prior_covariance;
      }
      double p2 = std::exp(-0.5*sum);
      return p1*p2;
   }
}



void RJMCMC::forward_modeling(std::vector<double> z, 
                        std::vector<double> sigma, 
                        struct PredictedResponse& predicted_response, 
                        double& misfit, double& rms)
{
   // Transform rjmcmc models (z in m, sigma in log10 scale) to radius in m, thickness in m, conductivity in S/m
   int n_layers = sigma.size()+1; // +1 means the Earth's core
   std::vector<double> radius(n_layers); // radius of top surface, in m
   std::vector<double> conductivity(n_layers); // in linear scale
   std::vector<double> thickness(n_layers); // thickness of each layer, in m
   std::vector<double> resistivity(n_layers);
   radius[0] = R0; // Earth's surface
   conductivity[0] = std::pow(10, sigma[0]);
   resistivity[0] = 1.0/conductivity[0];
   for (int i=1; i<n_layers-1; i++)
   {
      radius[i] = R0-z[i-1];
      conductivity[i] = std::pow(10, sigma[i]);
      thickness[i-1] = radius[i-1]-radius[i];
      resistivity[i] = 1.0/conductivity[i];
   }
   radius[n_layers-1] = R0-z_cmb; // core-mantle boundary
   conductivity[n_layers-1] = sigma_core; // core
   thickness[n_layers-2] = radius[n_layers-2]-radius[n_layers-1];
   thickness[n_layers-1] = radius[n_layers-1];
   resistivity[n_layers-1] = 1.0/conductivity[n_layers-1];

   // Compute misfit
   if (Qn_Ndata>0)
   {
      forward_Qn(radius, conductivity, Qn_SH_degree, Qn_periods, predicted_response.Qn_responses);
   }
   if (Cn_Ndata>0)
   {
      forward_Cn(radius, conductivity, Cn_SH_degree, Cn_periods, predicted_response.Cn_responses);
   }
   if (Dst_G2LTF_Ndata>0)
   {
      forward_G2LTF(radius, conductivity, Dst_G2LTF_SH_degree, Dst_G2LTF_SH_order, Dst_G2LTF_periods, 
                    GM_longitude, GM_colatitude, predicted_response.Dst_G2LTFs);
   }
   if (Sq_G2LTF_Ndata>0)
   {
      forward_G2LTF(radius, conductivity, Sq_G2LTF_SH_degree, Sq_G2LTF_SH_order, Sq_G2LTF_periods, 
                    GG_longitude, GG_colatitude, predicted_response.Sq_G2LTFs);
   }
   if (MT_appres_Ndata>0||MT_phase_Ndata>0)
   {
      std::vector<double> MT_all_appres;
      std::vector<double> MT_all_phase;
      forward_MT(thickness, resistivity, MT_all_periods, MT_all_appres, MT_all_phase);
      std::map<int,double> MT_all_periodid_appres;
      std::map<int,double> MT_all_periodid_phase;
      int idx = 0;
      for (std::map<int,double>::iterator i=MT_period_info.begin(); i!=MT_period_info.end(); i++)
      {
         int period_id = i->first;
         MT_all_periodid_appres.insert(std::make_pair(period_id,MT_all_appres[idx]));
         MT_all_periodid_phase.insert(std::make_pair(period_id,MT_all_phase[idx]));
         idx++;
      }
      // get real appres
      for (int i=0; i<MT_appres_Ndata; i++)
      {
         std::map<int, double>::iterator it = MT_all_periodid_appres.find(MT_appres_period_id[i]);
         assert(it!=MT_all_periodid_appres.end());
         predicted_response.MT_appres.push_back(it->second);
      }
      // get real phase
      for (int i=0; i<MT_phase_Ndata; i++)
      {
         std::map<int, double>::iterator it = MT_all_periodid_phase.find(MT_phase_period_id[i]);
         assert(it!=MT_all_periodid_phase.end());
         predicted_response.MT_phase.push_back(it->second);
      }
   }

   misfit = 0.;
   for (int i=0; i<Qn_Ndata; i++)
   {
      misfit += std::pow(std::abs(Qn_responses[i]-predicted_response.Qn_responses[i])/Qn_data_uncertainty[i], 2);
   }
   for (int i=0; i<Cn_Ndata; i++)
   {
      misfit += std::pow(std::abs(Cn_responses[i]-predicted_response.Cn_responses[i])/Cn_data_uncertainty[i], 2);
   }
   for (int i=0; i<Dst_G2LTF_Ndata; i++)
   {
      misfit += std::pow(std::abs(Dst_G2LTF_responses[i]-predicted_response.Dst_G2LTFs[i])/Dst_G2LTF_data_uncertainty[i], 2);
   }
   for (int i=0; i<Sq_G2LTF_Ndata; i++)
   {
      misfit += std::pow(std::abs(Sq_G2LTF_responses[i]-predicted_response.Sq_G2LTFs[i])/Sq_G2LTF_data_uncertainty[i], 2);
   }
   for (int i=0; i<MT_appres_Ndata; i++)
   {
      misfit += std::pow((MT_appres[i]-predicted_response.MT_appres[i])/MT_appres_uncertainty[i], 2);
   }
   for (int i=0; i<MT_phase_Ndata; i++)
   {
      misfit += std::pow((MT_phase[i]-predicted_response.MT_phase[i])/MT_phase_uncertainty[i], 2);
   }
   rms = std::sqrt(misfit/(MT_appres_Ndata+MT_phase_Ndata+Qn_Ndata+Cn_Ndata+Dst_G2LTF_Ndata+Sq_G2LTF_Ndata));
}



void RJMCMC::birth(std::vector<double> z, 
                   std::vector<double> sigma, 
                   std::vector<double>& z_new, 
                   std::vector<double>& sigma_new,
                   double& log_proposal_ratio,
                   bool& violate_prior)
{
   // Step 1: propose new interface
   int k = sigma.size();
   double z_proposal = z_min + (z_max-z_min)*uniform_distribution(random_engine);
   std::vector<double> z_vec(k+1,0.0);
   z_vec[0] = 0;
   z_vec[k] = mcmc_z_cmb;
   int j = 0;
   for (int i=1; i<k; i++)
   {
      z_vec[i] = z[j];
      j++;
   }
   int l = find_position(z_vec, z_proposal);

   // ensure that no layers thinner than h_min
   double z_up = z_vec[l];
   double z_low = z_vec[l+1];
   if (z_proposal-z_up<h_min || z_low-z_proposal<h_min)
   {
      violate_prior = true;
      return;
   }

   // update z
   z_new.resize(z.size()+1);
   for (int i=0; i<l; i++) z_new[i] = z[i];
   z_new[l] = z_proposal;
   for (int i=l+1; i<z_new.size(); i++) z_new[i] = z[i-1];


   // Step 2: propose new conductivity model
   // propose new conductivity
   double mu = sigma[l]; // mean model
   double gaussian_draw = gaussian_distribution(random_engine); // gaussian_draw = (sigma_proposal-mu)/std
   double sigma_proposal = mu + gaussian_draw*birth_death_std;
   // ensure that no conductivity is beyond minimum and maximum limits
   if (sigma_proposal<sigma_min || sigma_proposal>sigma_max)
   {
      violate_prior = true;
      return;
   }

   // chose upper or lower layer to be new layer
   if (uniform_distribution(random_engine) > 0.5) l++;

   // update conductivity
   sigma_new.resize(sigma.size()+1);
   for (int i=0; i<l; i++) sigma_new[i] = sigma[i];
   sigma_new[l] = sigma_proposal;
   for (int i=l+1; i<sigma_new.size(); i++) sigma_new[i] = sigma[i-1];


   // Step 3: cpmpute log proposal ratio
   log_proposal_ratio = 0.5*std::log(2*pi) + std::log(birth_death_std) + 0.5*std::pow(gaussian_draw,2);
}



void RJMCMC::death(std::vector<double> z, 
                   std::vector<double> sigma, 
                   std::vector<double>& z_new, 
                   std::vector<double>& sigma_new,
                   double& log_proposal_ratio,
                   bool& violate_prior)
{
   // Step 1: choose one interface
   // randomly choose one interface (k-1 layer interfaces with id from 0,1,...,k-2)
   int k = sigma.size();
   std::uniform_int_distribution<unsigned> randi(0,k-2);
   int l = randi(random_engine);

   // update z
   z_new.resize(z.size()-1);
   for (int i=0; i<l; i++) z_new[i] = z[i];
   for (int i=l; i<z_new.size(); i++) z_new[i] = z[i+1];


   // Step 2: propose new conductivity model
   // // using gaussian distribution (warning: proposal_ratio may need to be modified)
   // double mu = 0.5*(sigma[l]+sigma[l+1]); // mean model
   // double gaussian_draw = gaussian_distribution(random_engine);
   // double sigma_proposal = mu + gaussian_draw*birth_death_std;
   // // ensure that no conductivity is beyond minimum and maximum limits
   // if (sigma_proposal<sigma_min || sigma_proposal>sigma_max)
   // {
   //    violate_prior = true;
   //    return;
   // }

   // choose upper or lower layer conductivity as new value
   double sigma_proposal;
   if (uniform_distribution(random_engine) > 0.5)
   {
      sigma_proposal = sigma[l+1];
   }
   else
   {
      sigma_proposal = sigma[l];
   }
   double mu = 0.5*(sigma[l]+sigma[l+1]); // mean model
   double gaussian_draw = (sigma_proposal-mu)/birth_death_std;

   // update conductivity
   sigma_new.resize(sigma.size()-1);
   for (int i=0; i<l; i++) sigma_new[i] = sigma[i];
   sigma_new[l] = sigma_proposal;
   for (int i=l+1; i<sigma_new.size(); i++) sigma_new[i] = sigma[i+1];


   // Step 3: cpmpute log proposal ratio
   log_proposal_ratio = -0.5*std::log(2*pi) - std::log(birth_death_std) - 0.5*std::pow(gaussian_draw,2);
}



void RJMCMC::move(std::vector<double> z, 
                  std::vector<double> sigma, 
                  std::vector<double>& z_new, 
                  std::vector<double>& sigma_new,
                  bool& violate_prior)
{
   // randomly choose one interface (k-1 layers with id from 0,1,...,k-2)
   int k = sigma.size();
   std::uniform_int_distribution<unsigned> randi(0,k-2);
   int l = randi(random_engine);

   // propose a new depth
   // 1: move within [-hmin,hmin]
   // double z_proposal = z[l] + h_min*(uniform_distribution(random_engine)*2-1);

   // 2: 2021/11/12, move within [-move_std,move_std]
   // double z_proposal = z[l] + move_std*(uniform_distribution(random_engine)*2-1);

   // 3: 2022/06/15, using gaussian distribution to perturb depth
   double z_proposal = z[l] + gaussian_distribution(random_engine)*move_std;

   if (z_proposal<z_min || z_proposal>z_max)
   {
      violate_prior = true;
      return;
   }

   // find the new depth in which layer
   std::vector<double> z_vec(k+1,0.0);
   z_vec[0] = 0;
   z_vec[k] = mcmc_z_cmb;
   int j = 0;
   for (int i=1; i<k; i++)
   {
      z_vec[i] = z[j];
      j++;
   }
   int l_new = find_position(z_vec, z_proposal);
   if (l_new!=l && l_new!=l+1)
   {
      violate_prior = true;
      return;
   }

   // ensure that no layers thinner than h_min
   double z_up = z_vec[l];
   double z_low = z_vec[l+2];
   if (z_proposal-z_up<h_min || z_low-z_proposal<h_min)
   {
      violate_prior = true;
      return;
   }

   // update z
   z_new = z;
   z_new[l] = z_proposal;

   // update conductivity
   sigma_new = sigma;
}



void RJMCMC::update(std::vector<double> z, 
                    std::vector<double> sigma, 
                    std::vector<double>& z_new, 
                    std::vector<double>& sigma_new,
                    bool& violate_prior)
{
   // Step 1: update z
   z_new = z;
   sigma_new = sigma;

   // Step 2: propose new conductivity model
   // randomly choose one layer (k layers with id from 0,1,...,k-1)
   int k = sigma.size();
   std::uniform_int_distribution<unsigned> randi(0,k-1);
   int l = randi(random_engine);
   double sigma_proposal = sigma[l] + gaussian_distribution(random_engine)*update_std;

   // ensure that no conductivity is beyond minimum and maximum limits
   if (sigma_proposal<sigma_min || sigma_proposal>sigma_max)
   {
      violate_prior = true;
      return;
   }
   sigma_new[l] = sigma_proposal;
}



void RJMCMC::run_rjmcmc()
{
   double start_time,end_time;
   start_time = MPI_Wtime();

   // Samples
   std::vector<int> k_samples(maxit);
   std::vector< std::vector<double> > z_samples(maxit);
   std::vector< std::vector<double> > sigma_samples(maxit);
   std::vector< struct PredictedResponse > d_pred_samples(maxit);
   std::vector<double> rms_samples(maxit);

   // For acceptance rate
   int birth_count = 0, accept_birth_count = 0;
   int death_count = 0, accept_death_count = 0;
   int move_count = 0, accept_move_count = 0;
   int update_count = 0, accept_update_count = 0;

   // Initialization with a k_min layered model
   int k = k_min;
   std::vector<double> z(k-1);
   for (int i=0; i<k-1; i++)
   {
      z[i] = z_min+(z_max-z_min)*uniform_distribution(random_engine);
   }
   std::sort(z.begin(), z.end());
   std::vector<double> sigma(k);
   for (int i=0; i<k; i++)
   {
      sigma[i] = sigma_min+(sigma_max-sigma_min)*uniform_distribution(random_engine);
   }

   // Get initial misfit
   struct PredictedResponse predicted_response_old;
   double misfit_old, rms_old;
   forward_modeling(z, sigma, predicted_response_old, misfit_old, rms_old);
   double p_prior_old = get_prior_probability(z, sigma);

   // Probability of birth, death, move, update
   double prob[4] = {birth_prob, death_prob, move_prob, update_prob};
   std::vector<double> interval(5,0.0);
   double sum = 0.0;
   for (int i=0; i<4; i++)
   {
      sum += prob[i];
      interval[i+1] = sum;
   }

   // satrt rjmcmc iterations
   if (my_rank==0)
   {
      std::cout << "\n   ======================================== Iteration ========================================\n";
   }
   for (int iter=0; iter<maxit; iter++)
   {
      int temp_iter = iter+1;
      if (iter==0 && my_rank==0)
      {
         std::cout << "                       Iteration: " << std::setw(8) << iter 
                   << "      RMS misfit: " << std::setw(10) << rms_old << "\n";
      }
      if (temp_iter%print_every_iteration==0 && my_rank==0)
      {
         std::cout << "                       Iteration: " << std::setw(8) << temp_iter 
                   << "      RMS misfit: " << std::setw(10) << rms_old << "\n";
      }
      
      double rand_number = uniform_distribution(random_engine);
      int option = find_position(interval, rand_number);

      // Birth
      if (option==0)
      {
         birth_count++;

         // initialize the samples
         k_samples[iter] = k;
         z_samples[iter] = z;
         sigma_samples[iter] = sigma;
         d_pred_samples[iter] = predicted_response_old;
         rms_samples[iter] = rms_old;

         if (k<k_max)
         {
            std::vector<double> z_new;
            std::vector<double> sigma_new;
            double log_proposal_ratio;
            bool violate_prior = false;
            birth(z,sigma,z_new,sigma_new,log_proposal_ratio,violate_prior);

            if (!violate_prior)
            {
               // compute acceptance probability
               struct PredictedResponse predicted_response_new;
               double misfit_new, rms_new;
               forward_modeling(z_new, sigma_new, predicted_response_new, misfit_new, rms_new);
               double p_prior_new = get_prior_probability(z_new, sigma_new);
               double ln_alpha = std::log(p_prior_new)-std::log(p_prior_old)
                                 - 0.5*misfit_new + 0.5*misfit_old 
                                 + log_proposal_ratio;

               // accept or reject
               double r = std::log(uniform_distribution(random_engine));
               if (r < ln_alpha)
               {
                  accept_birth_count++;

                  k = k+1;
                  z = z_new;
                  sigma = sigma_new;
                  predicted_response_old = predicted_response_new;
                  misfit_old = misfit_new;
                  rms_old = rms_new;
                  p_prior_old = p_prior_new;

                  k_samples[iter] = k;
                  z_samples[iter] = z;
                  sigma_samples[iter] = sigma;
                  d_pred_samples[iter] = predicted_response_old;
                  rms_samples[iter] = rms_old;
               }
            }
         }
      }


      // Death
      else if (option==1)
      {
         death_count++;

         // initialize the samples
         k_samples[iter] = k;
         z_samples[iter] = z;
         sigma_samples[iter] = sigma;
         d_pred_samples[iter] = predicted_response_old;
         rms_samples[iter] = rms_old;

         if (k>k_min)
         {
            std::vector<double> z_new;
            std::vector<double> sigma_new;
            double log_proposal_ratio;
            bool violate_prior = false;
            death(z,sigma,z_new,sigma_new,log_proposal_ratio,violate_prior);

            if (!violate_prior)
            {
               // compute acceptance probability
               struct PredictedResponse predicted_response_new;
               double misfit_new, rms_new;
               forward_modeling(z_new, sigma_new, predicted_response_new, misfit_new, rms_new);
               double p_prior_new = get_prior_probability(z_new, sigma_new);
               double ln_alpha = std::log(p_prior_new)-std::log(p_prior_old)
                                 - 0.5*misfit_new + 0.5*misfit_old 
                                 + log_proposal_ratio;

               // accept or reject
               double r = std::log(uniform_distribution(random_engine));
               if (r < ln_alpha)
               {
                  accept_death_count++;

                  k = k-1;
                  z = z_new;
                  sigma = sigma_new;
                  predicted_response_old = predicted_response_new;
                  misfit_old = misfit_new;
                  rms_old = rms_new;
                  p_prior_old = p_prior_new;

                  k_samples[iter] = k;
                  z_samples[iter] = z;
                  sigma_samples[iter] = sigma;
                  d_pred_samples[iter] = predicted_response_old;
                  rms_samples[iter] = rms_old;
               }
            }
         }
      }


      // Move
      else if (option==2)
      {
         move_count++;

         // initialize the samples
         k_samples[iter] = k;
         z_samples[iter] = z;
         sigma_samples[iter] = sigma;
         d_pred_samples[iter] = predicted_response_old;
         rms_samples[iter] = rms_old;

         std::vector<double> z_new;
         std::vector<double> sigma_new;
         bool violate_prior = false;
         move(z,sigma,z_new,sigma_new,violate_prior);

         if (!violate_prior)
         {
            // compute acceptance probability
            struct PredictedResponse predicted_response_new;
            double misfit_new, rms_new;
            forward_modeling(z_new, sigma_new, predicted_response_new, misfit_new, rms_new);
            double p_prior_new = get_prior_probability(z_new, sigma_new);
            double ln_alpha = std::log(p_prior_new)-std::log(p_prior_old) - 0.5*misfit_new + 0.5*misfit_old;

            // accept or reject
            double r = std::log(uniform_distribution(random_engine));
            if (r < ln_alpha)
            {
               accept_move_count++;

               k = k;
               z = z_new;
               sigma = sigma_new;
               predicted_response_old = predicted_response_new;
               misfit_old = misfit_new;
               rms_old = rms_new;
               p_prior_old = p_prior_new;

               k_samples[iter] = k;
               z_samples[iter] = z;
               sigma_samples[iter] = sigma;
               d_pred_samples[iter] = predicted_response_old;
               rms_samples[iter] = rms_old;
            }
         }
      }


      // Update
      else if (option==3)
      {
         update_count++;

         // initialize the samples
         k_samples[iter] = k;
         z_samples[iter] = z;
         sigma_samples[iter] = sigma;
         d_pred_samples[iter] = predicted_response_old;
         rms_samples[iter] = rms_old;

         std::vector<double> z_new;
         std::vector<double> sigma_new;
         bool violate_prior = false;
         update(z,sigma,z_new,sigma_new,violate_prior);

         if (!violate_prior)
         {
            // compute acceptance probability
            struct PredictedResponse predicted_response_new;
            double misfit_new, rms_new;
            forward_modeling(z_new, sigma_new, predicted_response_new, misfit_new, rms_new);
            double p_prior_new = get_prior_probability(z_new, sigma_new);
            double ln_alpha = std::log(p_prior_new)-std::log(p_prior_old) - 0.5*misfit_new + 0.5*misfit_old;

            // accept or reject
            double r = std::log(uniform_distribution(random_engine));
            if (r < ln_alpha)
            {
               accept_update_count++;

               k = k;
               z = z_new;
               sigma = sigma_new;
               predicted_response_old = predicted_response_new;
               misfit_old = misfit_new;
               rms_old = rms_new;
               p_prior_old = p_prior_new;

               k_samples[iter] = k;
               z_samples[iter] = z;
               sigma_samples[iter] = sigma;
               d_pred_samples[iter] = predicted_response_old;
               rms_samples[iter] = rms_old;
            }
         }
      }
   }

   end_time = MPI_Wtime();
   if (my_rank==0) std::cout << "   Run time: " << end_time-start_time << " (s)\n";
   

   // Post-processing and saving results
   start_time = MPI_Wtime();
   if (my_rank==0) std::cout << "   Post-processing and saving results...\n";

   // Save acceptance rate
   std::ostringstream acceptance_rate_filename;
   acceptance_rate_filename << output_dir << "/acceptance_rate." 
                            << std::setfill('0') << std::setw(6) << my_rank << ".dat";
   std::ofstream acceptance_rate_stream(acceptance_rate_filename.str());
   acceptance_rate_stream << std::setw(7) << "Perturb" 
                           << std::setw(12) << "Count" 
                           << std::setw(25) << "Perturb probability" 
                           << std::setw(22) << "Acceptance count"
                           << std::setw(22) << "Acceptance rate\n";
   acceptance_rate_stream << std::setw(7) << "Birth" 
                           << std::setw(12) << birth_count
                           << std::setw(19) << (double)birth_count/(double)maxit
                           << std::setw(22) << accept_birth_count
                           << std::setw(22) << (double)accept_birth_count/(double)birth_count << "\n";
   acceptance_rate_stream << std::setw(7) << "Death" 
                           << std::setw(12) << death_count
                           << std::setw(19) << (double)death_count/(double)maxit
                           << std::setw(22) << accept_death_count
                           << std::setw(22) << (double)accept_death_count/(double)death_count << "\n";
   acceptance_rate_stream << std::setw(7) << "Move" 
                           << std::setw(12) << move_count
                           << std::setw(19) << (double)move_count/(double)maxit
                           << std::setw(22) << accept_move_count
                           << std::setw(22) << (double)accept_move_count/(double)move_count << "\n";
   acceptance_rate_stream << std::setw(7) << "Update" 
                           << std::setw(12) << update_count
                           << std::setw(19) << (double)update_count/(double)maxit
                           << std::setw(22) << accept_update_count
                           << std::setw(22) << (double)accept_update_count/(double)update_count << "\n";
   int accept_total_count = accept_birth_count+accept_death_count+accept_move_count+accept_update_count;
   acceptance_rate_stream << std::setw(7) << "Total" 
                           << std::setw(12) << maxit
                           << std::setw(19) << (double)1.0
                           << std::setw(22) << accept_total_count
                           << std::setw(22) << (double)accept_total_count/(double)maxit << "\n";
   acceptance_rate_stream.close();

   // Save predicted responses
   if (Qn_Ndata>0&&save_predicted_responses)
   {
      std::ostringstream out_filename;
      out_filename << output_dir << "/predicted_Qn_responses." 
                   << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream d_pred_stream(out_filename.str());
      for (int iter=burnin+1; iter<maxit; iter+=thin)
      {
         d_pred_stream << iter << "\t";
         for (int i=0; i<d_pred_samples[iter].Qn_responses.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Qn_responses[i].real() << "\t";
         }
         for (int i=0; i<d_pred_samples[iter].Qn_responses.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Qn_responses[i].imag() << "\t";
         }
         
         d_pred_stream << "\n";
      }
      d_pred_stream.close();
   }
   
   if (Cn_Ndata>0&&save_predicted_responses)
   {
      std::ostringstream out_filename;
      out_filename << output_dir << "/predicted_Cn_responses." 
                   << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream d_pred_stream(out_filename.str());
      for (int iter=burnin+1; iter<maxit; iter+=thin)
      {
         d_pred_stream << iter << "\t";
         for (int i=0; i<d_pred_samples[iter].Cn_responses.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Cn_responses[i].real() << "\t";
         }
         for (int i=0; i<d_pred_samples[iter].Cn_responses.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Cn_responses[i].imag() << "\t";
         }
         
         d_pred_stream << "\n";
      }
      d_pred_stream.close();
   }

   if (Dst_G2LTF_Ndata>0&&save_predicted_responses)
   {
      std::ostringstream out_filename;
      out_filename << output_dir << "/predicted_Dst_G2LTF_responses." 
                   << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream d_pred_stream(out_filename.str());
      for (int iter=burnin+1; iter<maxit; iter+=thin)
      {
         d_pred_stream << iter << "\t";
         for (int i=0; i<d_pred_samples[iter].Dst_G2LTFs.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Dst_G2LTFs[i].real() << "\t";
         }
         for (int i=0; i<d_pred_samples[iter].Dst_G2LTFs.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Dst_G2LTFs[i].imag() << "\t";
         }
         
         d_pred_stream << "\n";
      }
      d_pred_stream.close();
   }

   if (Sq_G2LTF_Ndata>0&&save_predicted_responses)
   {
      std::ostringstream out_filename;
      out_filename << output_dir << "/predicted_Sq_G2LTF_responses." 
                   << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream d_pred_stream(out_filename.str());
      for (int iter=burnin+1; iter<maxit; iter+=thin)
      {
         d_pred_stream << iter << "\t";
         for (int i=0; i<d_pred_samples[iter].Sq_G2LTFs.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Sq_G2LTFs[i].real() << "\t";
         }
         for (int i=0; i<d_pred_samples[iter].Sq_G2LTFs.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].Sq_G2LTFs[i].imag() << "\t";
         }
         
         d_pred_stream << "\n";
      }
      d_pred_stream.close();
   }

   if (MT_appres_Ndata>0&&save_predicted_responses)
   {
      std::ostringstream out_filename;
      out_filename << output_dir << "/predicted_MT_appres_responses." 
                   << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream d_pred_stream(out_filename.str());
      for (int iter=burnin+1; iter<maxit; iter+=thin)
      {
         d_pred_stream << iter << "\t";
         for (int i=0; i<d_pred_samples[iter].MT_appres.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].MT_appres[i] << "\t";
         }
         d_pred_stream << "\n";
      }
      d_pred_stream.close();
   }

   if (MT_phase_Ndata>0&&save_predicted_responses)
   {
      std::ostringstream out_filename;
      out_filename << output_dir << "/predicted_MT_phase_responses." 
                   << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream d_pred_stream(out_filename.str());
      for (int iter=burnin+1; iter<maxit; iter+=thin)
      {
         d_pred_stream << iter << "\t";
         for (int i=0; i<d_pred_samples[iter].MT_phase.size(); i++)
         {
            d_pred_stream << d_pred_samples[iter].MT_phase[i] << "\t";
         }
         d_pred_stream << "\n";
      }
      d_pred_stream.close();
   }

   // Save RMS
   if (save_rms_misfit)
   {
      std::ostringstream rms_filename;
      rms_filename << output_dir << "/rms." << std::setfill('0') << std::setw(6) << my_rank << ".dat";
      std::ofstream rms_stream(rms_filename.str());
      for (int iter=0; iter<maxit; iter++)
      {
         rms_stream << iter << "\t" << rms_samples[iter] << "\n";
      }
      rms_stream.close();
   }

   /* Parallel post-processing */
   // Keep one in every thin samples
   std::vector<int> keep_k_samples;
   std::vector<double> keep_z_samples;
   std::vector<double> keep_sigma_samples;
   for (int iter=burnin+1; iter<maxit; iter+=thin)
   {
      keep_k_samples.push_back(k_samples[iter]);
      for (int i=0; i<z_samples[iter].size(); i++)
      {
         keep_z_samples.push_back(z_samples[iter][i]);
      }
      for (int i=0; i<sigma_samples[iter].size(); i++)
      {
         keep_sigma_samples.push_back(sigma_samples[iter][i]);
      }
   }
   // if (my_rank==0) std::cout << keep_k_samples.size() << "\t" << keep_z_samples.size() << "\t" << keep_sigma_samples.size() << "\n";

   // Send all samples to master processor
   int comm_size;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   int z_samples_num;
   int sigma_samples_num;
   int *msg_k_samples;
   double *msg_z_samples;
   double *msg_sigma_samples;
   if (my_rank==0)
   {
      std::vector<int> keep_k_samples_all;
      std::vector< std::vector<double> > keep_z_samples_all;
      std::vector< std::vector<double> > keep_sigma_samples_all;
      int z_idx = 0, sigma_idx = 0;
      for (int s=0; s<keep_k_samples.size(); s++)
      {
         keep_k_samples_all.push_back(keep_k_samples[s]);

         std::vector<double> z_sample;
         for (int i=0; i<keep_k_samples[s]-1; i++)
         {
            z_sample.push_back(keep_z_samples[z_idx]);
            z_idx++;
         }
         keep_z_samples_all.push_back(z_sample);

         std::vector<double> sigma_sample;
         for (int i=0; i<keep_k_samples[s]; i++)
         {
            sigma_sample.push_back(keep_sigma_samples[sigma_idx]);
            sigma_idx++;
         }
         keep_sigma_samples_all.push_back(sigma_sample);
      }
      // std::cout << z_idx << "\t" << sigma_idx << "\n";

      for (int p=1; p<comm_size; p++)
      {
         // receive k
         msg_k_samples = new int[keep_k_samples.size()];
         MPI_Recv(msg_k_samples, keep_k_samples.size(), MPI_INT, p, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         
         // receive z
         MPI_Recv(&z_samples_num, 1, MPI_INT, p, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         msg_z_samples = new double[z_samples_num];
         MPI_Recv(msg_z_samples, z_samples_num, MPI_DOUBLE, p, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         // receive sigma
         MPI_Recv(&sigma_samples_num, 1, MPI_INT, p, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         msg_sigma_samples = new double[sigma_samples_num];
         MPI_Recv(msg_sigma_samples, sigma_samples_num, MPI_DOUBLE, p, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         z_idx = 0, sigma_idx = 0;
         for (int s=0; s<keep_k_samples.size(); s++)
         {
            keep_k_samples_all.push_back(msg_k_samples[s]);

            std::vector<double> z_sample;
            for (int i=0; i<msg_k_samples[s]-1; i++)
            {
               z_sample.push_back(msg_z_samples[z_idx]);
               z_idx++;
            }
            keep_z_samples_all.push_back(z_sample);

            std::vector<double> sigma_sample;
            for (int i=0; i<msg_k_samples[s]; i++)
            {
               sigma_sample.push_back(msg_sigma_samples[sigma_idx]);
               sigma_idx++;
            }
            keep_sigma_samples_all.push_back(sigma_sample);
         }

         delete msg_k_samples;
         delete msg_z_samples;
         delete msg_sigma_samples;
      }
      // std::cout << keep_k_samples.size() << "\t" << keep_k_samples_all.size() << "\n";

      // Save results
      save_results(keep_k_samples_all, keep_z_samples_all, keep_sigma_samples_all);
   }
   else
   {
      // send k
      msg_k_samples = keep_k_samples.data();
      MPI_Send(msg_k_samples, keep_k_samples.size(), MPI_INT, 0, 101, MPI_COMM_WORLD);

      // send z
      z_samples_num = keep_z_samples.size();
      MPI_Send(&z_samples_num, 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
      msg_z_samples = keep_z_samples.data();
      MPI_Send(msg_z_samples, z_samples_num, MPI_DOUBLE, 0, 103, MPI_COMM_WORLD);

      // send sigma
      sigma_samples_num = keep_sigma_samples.size();
      MPI_Send(&sigma_samples_num, 1, MPI_INT, 0, 104, MPI_COMM_WORLD);
      msg_sigma_samples = keep_sigma_samples.data();
      MPI_Send(msg_sigma_samples, sigma_samples_num, MPI_DOUBLE, 0, 105, MPI_COMM_WORLD);
   }

   end_time = MPI_Wtime();
   if (my_rank==0) std::cout << "   Run time: " << end_time-start_time << " (s)\n";
}



void RJMCMC::save_results(std::vector<int> k_samples,
                          std::vector< std::vector<double> > z_samples,
                          std::vector< std::vector<double> > sigma_samples)
{
   if (write_model_samples)
   {
      std::string model_samples_filename = output_dir + "/model_samples.dat";
      std::ofstream model_samples_stream(model_samples_filename);
      for (int i=0; i<k_samples.size(); i++)
      {
         model_samples_stream << k_samples[i] << "\t" << 0 << "\t";
         for (int j=0; j<z_samples[i].size(); j++)
         {
            model_samples_stream << z_samples[i][j] << "\t";
         }
         for (int j=0; j<sigma_samples[i].size(); j++)
         {
            model_samples_stream << sigma_samples[i][j] << "\t";
         }
         model_samples_stream << "\n";
      }
      model_samples_stream.close();
   }

   // Linspace of conductivity
   std::vector<double> sigma_linspace(n_sigma_blocks+1);
   double sigma_interval = delta_sigma/n_sigma_blocks;
   for (int i=0; i<sigma_linspace.size(); i++)
   {
      sigma_linspace[i] = sigma_min + i*sigma_interval;
   }
   std::string sigma_blocks_filename = output_dir + "/sigma_blocks.dat";
   std::ofstream sigma_blocks_stream(sigma_blocks_filename);
   for (int i=0; i<sigma_linspace.size(); i++)
   {
      sigma_blocks_stream << sigma_linspace[i] << "\n";
   }
   sigma_blocks_stream.close();
   

   // Linspace of depth
   std::vector<double> depth_linspace(n_depth_blocks+1);
   double depth_interval = (z_max-z_min)/n_depth_blocks;
   for (int i=0; i<depth_linspace.size(); i++)
   {
      depth_linspace[i] = z_min + i*depth_interval;
   }
   std::vector<double> lower_depths(n_depth_blocks);
   for (int i=0; i<lower_depths.size(); i++)
   {
      lower_depths[i] = depth_linspace[i+1];
   }
   std::string depth_blocks_filename = output_dir + "/depth_blocks.dat";
   std::ofstream depth_blocks_stream(depth_blocks_filename);
   for (int i=0; i<depth_linspace.size(); i++)
   {
      depth_blocks_stream << depth_linspace[i] << "\n";
   }
   depth_blocks_stream.close();

   // Conductivity values of each depth interval
   std::vector< std::vector<double> > sigma_vector_of_each_depth(n_depth_blocks);
   std::vector<int> interface_depth_count(n_depth_blocks+1,0);
   for (int s=0; s<k_samples.size(); s++)
   {
      int k = k_samples[s];
      std::vector<double> z_sample(k);
      std::vector<double> sigma_sample = sigma_samples[s];
      for (int i=0; i<k-1; i++)
      {
         z_sample[i] = z_samples[s][i];
      }
      z_sample[k-1] = mcmc_z_cmb+10; // +10 is more safe when comparing z_sample and lower_depths

      int j_start = 0;
      for (int i=0; i<z_sample.size(); i++)
      {
         double z = z_sample[i];
         double sigma = sigma_sample[i];
         for (int j=j_start; j<lower_depths.size(); j++) // start from bottom layer interface
         {
            if (lower_depths[j] <= z)
            {
               sigma_vector_of_each_depth[j].push_back(sigma);
            }
            else
            {
               interface_depth_count[j] += 1;
               j_start = j;
               break;
            }
         }
      }
   }

   // conductivity count matrix, n_depth_blocks x n_sigma_blocks
   std::vector< std::vector<int> > count_matrix(n_depth_blocks);
   std::vector<double> median_model(n_depth_blocks);
   std::vector<double> mean_model(n_depth_blocks);
   for (int i=0; i<n_depth_blocks; i++)
   {
      // mean and median models
      std::vector<double> sigma_vector = sigma_vector_of_each_depth[i];
      std::sort(sigma_vector.begin(), sigma_vector.end());
      int size = sigma_vector.size();
      mean_model[i] = std::accumulate(sigma_vector.begin(), sigma_vector.end(), 0.0)/(double)size;
      if (size%2 == 0)
      {
         int loc = size/2;
         median_model[i] = (sigma_vector[loc]+sigma_vector[loc-1])/2.0;
      }
      else
      {
         int loc = (size-1)/2;
         median_model[i] = sigma_vector[loc];
      }

      // conductivity count matrix
      count_matrix[i].resize(n_sigma_blocks);
      for (int j=0; j<sigma_vector.size(); j++)
      {
         int pos = find_position(sigma_linspace, sigma_vector[j]);
         // assert(pos>=0 && pos<n_sigma_blocks);
         count_matrix[i][pos] += 1;
      }
   }

   // Probability density function matrix (conductivity distribution), n_depth_blocks x n_sigma_blocks
   std::vector< std::vector<double> > pdf_matrix(n_depth_blocks);
   std::vector< std::vector<double> > pdf_matrix_normalize(n_depth_blocks);
   std::vector<double> most_probable_model(n_depth_blocks);
   std::vector<double> left_credible_model(n_depth_blocks);
   std::vector<double> right_credible_model(n_depth_blocks);
   for (int i=0; i<n_depth_blocks; i++)
   {
      std::vector<int> count_matrix_row = count_matrix[i];
      int total = std::accumulate(count_matrix_row.begin(), count_matrix_row.end(), 0);
      pdf_matrix[i].resize(n_sigma_blocks);
      for (int j=0; j<n_sigma_blocks; j++)
      {
         pdf_matrix[i][j] =  (double)count_matrix[i][j]/(double)total;
      }

      auto max_pdf_pos = std::max_element(pdf_matrix[i].begin(), pdf_matrix[i].end());
      pdf_matrix_normalize[i].resize(n_sigma_blocks);
      for (int j=0; j<n_sigma_blocks; j++)
      {
         pdf_matrix_normalize[i][j] = pdf_matrix[i][j]/(*max_pdf_pos);
      }
      int most_probable_pos =  max_pdf_pos-pdf_matrix[i].begin();
      most_probable_model[i] = sigma_linspace[most_probable_pos]; // using left one

      double prob_sum = 0.0;
      int left_pos = 0; // sigma_linspace[0] in default
      int right_pos = sigma_linspace.size()-1; // sigma_linspace[sigma_linspace.size()-1] in default
      for (int j=0; j<n_sigma_blocks; j++)
      {
         prob_sum += pdf_matrix[i][j];
         if (prob_sum < left_credible_interval)
         {
            left_pos = j; // keep updating until prob_sum >= left_credible_interval
         }
         else if (prob_sum > credible_interval)
         {
            right_pos = j;
            break;
         }
      }
      left_credible_model[i] = sigma_linspace[left_pos];
      right_credible_model[i] = sigma_linspace[right_pos];
   }
   std::string mean_out_filename = output_dir + "/mean_model.dat";
   std::ofstream mean_stream(mean_out_filename);
   std::string median_out_filename = output_dir + "/median_model.dat";
   std::ofstream median_stream(median_out_filename);
   std::string most_probable_out_filename = output_dir + "/most_probable_model.dat";
   std::ofstream most_probable_stream(most_probable_out_filename);
   std::string credible_out_filename = output_dir + "/credible_interval.dat";
   std::ofstream credible_stream(credible_out_filename);
   std::string pdf_matrix_filename = output_dir + "/pdf_matrix.dat";
   std::ofstream pdf_matrix_stream(pdf_matrix_filename);
   std::string pdf_matrix_normalize_filename = output_dir + "/pdf_matrix_normalize.dat";
   std::ofstream pdf_matrix_normalize_stream(pdf_matrix_normalize_filename);
   mean_stream.setf(std::ios::scientific);
   median_stream.setf(std::ios::scientific);
   most_probable_stream.setf(std::ios::scientific);
   credible_stream.setf(std::ios::scientific);
   for (int i=0; i<n_depth_blocks; i++)
   {
      mean_stream << std::setw(16)  << depth_linspace[i] << std::setw(16)  << mean_model[i] << "\n";
      median_stream << std::setw(16)  << depth_linspace[i] << std::setw(16)  << median_model[i] << "\n";
      most_probable_stream << std::setw(16)  << depth_linspace[i] << std::setw(16)  << most_probable_model[i] << "\n";
      credible_stream << std::setw(16)  << depth_linspace[i] << std::setw(16)  << left_credible_model[i] << std::setw(16)  << right_credible_model[i] << "\n";
      for (int j=0; j<n_sigma_blocks; j++)
      {
         pdf_matrix_stream << pdf_matrix[i][j] << "\t";
         pdf_matrix_normalize_stream << pdf_matrix_normalize[i][j] << "\t";
      }
      pdf_matrix_stream << "\n";
      pdf_matrix_normalize_stream << "\n";
   }
   mean_stream.close();
   median_stream.close();
   most_probable_stream.close();
   credible_stream.close();
   pdf_matrix_stream.close();
   pdf_matrix_normalize_stream.close();
   
   // Interface depth probability 
   std::vector<double> interface_depths(n_depth_blocks);
   for (int i=0; i<n_depth_blocks; i++)
   {
      interface_depths[i] = 0.5*(depth_linspace[i]+depth_linspace[i+1]);
   }
   std::vector<double> interface_depth_probability(n_depth_blocks);
   int total_size = std::accumulate(interface_depth_count.begin(), interface_depth_count.end(), 0);
   for (int i=0; i<n_depth_blocks; i++)
   {
      interface_depth_probability[i] = (double)interface_depth_count[i]/(double)total_size;
   }
   std::string interface_depth_prob_filename = output_dir + "/interface_depth_probability.dat";
   std::ofstream interface_depth_prob_stream(interface_depth_prob_filename);
   for (int i=0; i<interface_depth_probability.size(); i++)
   {
      interface_depth_prob_stream << interface_depths[i] << "\t" << interface_depth_probability[i] << "\n";
   }
   interface_depth_prob_stream.close();

   // Number of layer probability
   std::vector<int> layer_id(k_max-k_min+1);
   std::vector<int> num_of_layer_count(k_max-k_min+1);
   std::vector<double> num_of_layer_probability(k_max-k_min+1);
   for (int i=0; i<layer_id.size(); i++)
   {
      layer_id[i] = k_min + i;
      num_of_layer_count[i] = std::count(k_samples.begin(), k_samples.end(), layer_id[i]);
   }
   int total_layer_count = std::accumulate(num_of_layer_count.begin(), num_of_layer_count.end(), 0);
   for (int i=0; i<layer_id.size(); i++)
   {
      num_of_layer_probability[i] = (double)num_of_layer_count[i]/(double)total_layer_count;
   }
   std::string num_of_layer_prob_filename = output_dir + "/num_of_layer_probability.dat";
   std::ofstream num_of_layer_prob_stream(num_of_layer_prob_filename);
   for (int i=0; i<num_of_layer_probability.size(); i++)
   {
      num_of_layer_prob_stream << layer_id[i] << "\t" << num_of_layer_probability[i] << "\n";
   }
   num_of_layer_prob_stream.close();



   /* Get MCMC final median model with coarse depth */
   // Linspace of depth for final model
   std::vector<double> depth_linspace_coarse;
   if (model_interface_depth_file!="")
   {
      std::ifstream model_in_stream(model_interface_depth_file);
      std::vector<double> depth;
      std::string line;
      while (getline(model_in_stream, line))
      {
         if (*(line.begin()) == '#' || line == "") continue;
         std::stringstream sstr(line);
         double temp_depth;
         sstr >> temp_depth;
         depth_linspace_coarse.push_back(temp_depth);
      }
      model_in_stream.close();
      n_depth_blocks_coarse = depth_linspace_coarse.size()-1;
   }
   else
   {
      depth_linspace_coarse.resize(n_depth_blocks_coarse+1);
      double depth_interval_coarse = (z_max-z_min)/n_depth_blocks_coarse;
      for (int i=0; i<depth_linspace_coarse.size(); i++)
      {
         depth_linspace_coarse[i] = z_min + i*depth_interval_coarse;
      }
   }
   
   std::vector<double> lower_depths_coarse(n_depth_blocks_coarse);
   for (int i=0; i<lower_depths_coarse.size(); i++)
   {
      lower_depths_coarse[i] = depth_linspace_coarse[i+1];
   }

   // Conductivity values of each depth interval
   std::vector< std::vector<double> > sigma_vector_of_each_depth_coarse(n_depth_blocks_coarse);
   for (int s=0; s<k_samples.size(); s++)
   {
      int k = k_samples[s];
      std::vector<double> z_sample(k);
      std::vector<double> sigma_sample = sigma_samples[s];
      for (int i=0; i<k-1; i++)
      {
         z_sample[i] = z_samples[s][i];
      }
      z_sample[k-1] = mcmc_z_cmb+10; // +10 is more safe when comparing z_sample and lower_depths

      int j_start = 0;
      for (int i=0; i<z_sample.size(); i++)
      {
         double z = z_sample[i];
         double sigma = sigma_sample[i];
         for (int j=j_start; j<lower_depths_coarse.size(); j++) // start from bottom layer interface
         {
            if (lower_depths_coarse[j] <= z)
            {
               sigma_vector_of_each_depth_coarse[j].push_back(sigma);
            }
            else
            {
               j_start = j;
               break;
            }
         }
      }
   }

   // Median model with lower depth resolution
   std::vector<double> median_model_final(n_depth_blocks_coarse); 
   for (int i=0; i<n_depth_blocks_coarse; i++)
   {
      std::vector<double> sigma_vector = sigma_vector_of_each_depth_coarse[i];
      std::sort(sigma_vector.begin(), sigma_vector.end());
      int size = sigma_vector.size();
      if (size%2 == 0)
      {
         int loc = size/2;
         median_model_final[i] = (sigma_vector[loc]+sigma_vector[loc-1])/2.0;
      }
      else
      {
         int loc = (size-1)/2;
         median_model_final[i] = sigma_vector[loc];
      }
   }

   // Save final model
   std::string median_stream_coarse_filename = output_dir + "/median_model_final.dat";
   std::ofstream median_stream_final(median_stream_coarse_filename);
   median_stream_final.setf(std::ios::scientific);
   for (int i=0; i<n_depth_blocks_coarse; i++)
   {
      median_stream_final << std::setw(16) << depth_linspace_coarse[i] 
                          << std::setw(16) << std::pow(10,median_model_final[i]) << "\n";
   }
   median_stream_final << std::setw(16) << z_cmb << std::setw(16) << sigma_core << "\n";
   median_stream_final.close();
}