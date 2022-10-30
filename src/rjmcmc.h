// Trans-dimensional Bayesian inversion of magnetotelluric and 
// geomagnetic depth sounding data

// Copyright (c) 2021-2022 
// Authors: Hongbo Yao and Zhengyong Ren
// Institute: Central South University (CSU)
// Email: yaohongbo@csu.edu.cn and renzhengyong@csu.edu.cn
// GitHub Page: https://github.com/hongbo-yao
// Researchgate Page: https://www.researchgate.net/profile/Hongbo_Yao2

#ifndef _RJMCMC_H
#define _RJMCMC_H

#include <vector>
#include <map>
#include <string>
#include "utils.h"

// Predicted responses
struct PredictedResponse
{
   std::vector<double>     MT_appres;
   std::vector<double>     MT_phase;
   std::vector<Dcomplex>   Qn_responses;
   std::vector<Dcomplex>   Cn_responses;
   std::vector<Dcomplex>   Dst_G2LTFs;
   std::vector<Dcomplex>   Sq_G2LTFs;
};


class RJMCMC
{
public:
   int my_rank;

   /* Observed data */
   std::string             observed_data_file = "";
   std::string             station_name = "";
   double                  GG_longitude = 9999;
   double                  GG_colatitude = 9999;
   double                  GM_longitude = 9999;
   double                  GM_colatitude = 9999;
   int                     N_data;

   // Magnetotelluric response function
   // MT data
   std::map<int,double>    MT_period_info;
   std::vector<double>     MT_all_periods;
   int                     MT_appres_Ndata = 0;
   std::vector<int>        MT_appres_period_id;
   std::vector<double>     MT_appres_periods;
   std::vector<double>     MT_appres;
   std::vector<double>     MT_appres_uncertainty;
   int                     MT_phase_Ndata = 0;
   std::vector<int>        MT_phase_period_id;
   std::vector<double>     MT_phase_periods;
   std::vector<double>     MT_phase;
   std::vector<double>     MT_phase_uncertainty;

   // Global GDS transfer function: Qn-response
   // Qn data (Dst+Sq)
   int                     Qn_Ndata = 0;
   std::vector<double>     Qn_periods;
   std::vector<double>     Qn_SH_degree;
   std::vector<Dcomplex>   Qn_responses;
   std::vector<double>     Qn_data_uncertainty;

   // Local GDS transfer function: Cn-response
   // Cn data (Dst+Sq)
   int                     Cn_Ndata = 0;
   std::vector<double>     Cn_periods;
   std::vector<double>     Cn_SH_degree;
   std::vector<Dcomplex>   Cn_responses;
   std::vector<double>     Cn_data_uncertainty;

   // Global-to-local (G2L) GDS transfer function
   // Dst G2L TFs
   int                     Dst_G2LTF_Ndata = 0;
   std::vector<double>     Dst_G2LTF_periods;
   std::vector<double>     Dst_G2LTF_SH_degree;
   std::vector<double>     Dst_G2LTF_SH_order;
   std::vector<Dcomplex>   Dst_G2LTF_responses;
   std::vector<double>     Dst_G2LTF_data_uncertainty;

   // Sq G2L TFs
   int                     Sq_G2LTF_Ndata = 0;
   std::vector<double>     Sq_G2LTF_periods;
   std::vector<double>     Sq_G2LTF_SH_degree;
   std::vector<double>     Sq_G2LTF_SH_order;
   std::vector<Dcomplex>   Sq_G2LTF_responses;
   std::vector<double>     Sq_G2LTF_data_uncertainty;


   /* Common RJMCMC parameters */
   // Minimum and Maximum number of layers
   int      k_min = 2, k_max = 40;
   // Depth of core-mantle boundary
   double   mcmc_z_cmb = z_cmb;
   // Minimum and Maximum depths to each interface
   double   z_min = 0, z_max = z_cmb;
   // Minimum and Maximum log10 conductivity
   // delta_sigma=sigma_max-sigma_min
   double   sigma_min = -4, sigma_max = 2, delta_sigma; 
   // Probability of bitrh, death, move, update
   double   birth_prob = 0.25, death_prob = 0.25, move_prob = 0.25, update_prob = 0.25;
   // Number of mcmc iterations
   int      maxit = 100000;
   // Print every iterations
   int      print_every_iteration = 1000;
   

   /* Priori parameters */
   double   bestfit_halfspace_model;
   // Minimun allowed layer thickness
   double   h_min;
   // Priori distribution type
   // 0: uniform distribution
   // 1: normal distribution with reference model constraint
   // 2: normal distribution with smooth constraint
   int      prior_type = 2;
   // Standard derivation for the reference model perturbation in log10 conductivity
   // works for prior_type=1
   double   log10_conductivity_reference_std = 2.0;
   // Standard derivation for the difference in layer log10 conductivity
   // works for prior_type=2
   double   log10_conductivity_gradient_std = 1.0;


   /* Proposal parameters */
   // Standard deviation of log10 conductivity for birth/death pertubation
   double   birth_death_std;
   // Standard deviation of log10 conductivity for update pertubation
   double   update_std;
   // Standard deviation of depth in m for move pertubation
   double   move_std;


   /* Post-processing */
   // Burn in period
   int      burnin = 10000;
   // Keep one samples every thin_samples
   int      thin = 10;
   // Number of blocks in depth, such as 200
   int      n_depth_blocks = 200;
   // Number of blocks in conductivity, such as 100
   int      n_sigma_blocks = 100;
   // Bayesian credibl interval
   double   credible_interval = 0.95; // such as 0.95
   double   left_credible_interval; // 1-credible_interval
   // Model interface depth file used to generate the final model
   std::string model_interface_depth_file = "";
   // Number of blocks in depth for final model, less than n_depth_blocks, such as 40, 
   // only used when no model_interface_depth_file input
   int      n_depth_blocks_coarse = 40;
   // Flag of writing predicted responses
   // 0: do not write responses; 1: only write responses of rank 0; 2: write responses of all processes, default=1
   int      write_predicted_responses = 1;
   // Flag of writing rms misfit
   // 0: do not write rms; 1: only write rms of rank 0; 2: write rms of all processes, default=1
   int      write_rms_misfit = 1;
   // Output directory, default = "output"
   std::string output_dir = "output";
   // Flag of writing model samples, 0: do not write; 1: write; default=0
   int      write_model_samples = 0;

   /* For current process */
   bool     save_predicted_responses = false;
   bool     save_rms_misfit = false;

public:
   RJMCMC(std::string config_file);
   ~RJMCMC();

public:
   // Set default parameters
   void initialize();

   // Read and setup parameters
   void read_input_parameters(std::string config_file);
   void print_computation_parameters(std::string out="Bayesian_parameters.log");

   // Read observed data
   void read_observed_data();

   // Compute data misfit and predicted response
   void forward_modeling(std::vector<double> z, 
                         std::vector<double> sigma, 
                         struct PredictedResponse& predicted_response, 
                         double& misfit, double& rms);

   // Compute best-fit halfspace model as reference/initial model
   double get_bestfit_halfspace_model();
   
   // Compute prior probability on conductivity
   double get_prior_probability(std::vector<double> z, std::vector<double> sigma);

   // Randomly add a new layer
   void birth(std::vector<double> z, 
              std::vector<double> sigma, 
              std::vector<double>& z_new, 
              std::vector<double>& sigma_new,
              double& log_proposal_ratio,
              bool& violate_prior);

   // Randomly delete a layer
   void death(std::vector<double> z, 
              std::vector<double> sigma, 
              std::vector<double>& z_new, 
              std::vector<double>& sigma_new,
              double& log_proposal_ratio,
              bool& violate_prior);

   // Randomly move one interface
   void move(std::vector<double> z, 
             std::vector<double> sigma, 
             std::vector<double>& z_new, 
             std::vector<double>& sigma_new,
             bool& violate_prior);

   // Update conductivity values of all layers
   void update(std::vector<double> z, 
               std::vector<double> sigma, 
               std::vector<double>& z_new, 
               std::vector<double>& sigma_new,
               bool& violate_prior);

   // Run RJMCMC algorithm to do trans-dimensional Bayesian inversion
   void run_rjmcmc();

   // Post-process MCMC samples and output results for plotting in Python or Matlab
   void save_results(std::vector<int> k_samples,
                     std::vector< std::vector<double> > z_samples,
                     std::vector< std::vector<double> > sigma_samples);
};

#endif // _RJMCMC_H 