
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_tm2.c -- Very basic implementation of the TMM
*
******************************************************************************/

#include "imd.h"
#include <complex.h>


#define TMM_t0_suggest //wenn ja, wird t0 so angepasst dass laser E-field gerade @ threshold für laser aktivieren

#ifdef TTM1D
#define node  l1[i]
#define node2 l2[i]
#else
#define node  l1[i][j][k]
#define node2 l2[i][j][k]
#endif

#define RHOMIN  30  //kg/m^3 (active cell abh. von minatoms und dens)
                    // NOTE: This might be too small but should be fine. 
                    // No extrapolation needed if this value is the same as the lowest density in the files

#ifdef TFM
#define tm2_get_epsilon tm2_get_epsilon_TFM
#endif
#ifdef TTT
#ifdef MATERIALAL // Al
#define tm2_get_epsilon tm2_get_epsilon_Al
#define calc_y_bar calc_y_bar_Al
#define calc_beta calc_beta_Al
#define calc_w calc_w_Al
#endif // MATERIALAL
#ifdef MATERIALCU // Cu
#define tm2_get_epsilon tm2_get_epsilon_Cu
#define calc_y_bar calc_y_bar_Cu
#define calc_beta calc_beta_Cu
#define calc_w calc_w_Cu
#endif // MATERIALCU
#endif


int tm2_init()
{
  alloc1darr(real,tm2_eps_real_arr_global,global_fd_dim.x);
  alloc1darr(real,tm2_eps_imag_arr_global,global_fd_dim.x);
  alloc1darr(real,tm2_alpha_arr_global,global_fd_dim.x);
  alloc1darr(real,tm2_alpha_int_arr_global,global_fd_dim.x);

  alloc1darr(real,tm2_eps_real_arr_local,global_fd_dim.x);
  alloc1darr(real,tm2_eps_imag_arr_local,global_fd_dim.x);
  alloc1darr(real,tm2_alpha_arr_local,global_fd_dim.x);
  
#ifdef TFM
  tm2_read_arr(eps_filename,&eps_bb_data,3,&eps_bb_rows);
  tm2_read_arr(k12_filename,&K12,2,&K12_rows); //K1 und K1 num.eval.integrale
#endif

  alloc1darr(int,tm2_active_cell_global,global_fd_dim.x);
  alloc1darr(int,tm2_active_cell_local,global_fd_dim.x);

  alloc1darr(double,tm2_Qabs,global_fd_dim.x);
  alloc1darr(double,tm2_Qabs_scat,global_fd_dim.x);
  int i;

  tm2_dt=timestep*10.18/1e15; // in sek.
  for(i=0;i<global_fd_dim.x;++i)
  {
    tm2_eps_real_arr_local[i]=tm2_eps_real_arr_global[i]=0.0;
    tm2_eps_imag_arr_local[i]=tm2_eps_imag_arr_global[i]=0.0;
    tm2_alpha_arr_local[i]=tm2_alpha_arr_global[i]=0.0;
    tm2_active_cell_local[i]=tm2_active_cell_global[i]=0;

    tm2_Qabs[i]=0.0;
    tm2_Qabs_scat[i]=0.0;
  }

  real freq=c0/lambda;
  omega_laser=2.0*M_PI*freq;
  k0=2.0*M_PI/lambda;
  k0/=1e10; //weil fd_h.x in Angstrom

#ifdef TMM_t0_suggest
  // verschiebe t0 sodasss laser E-feld gerade threshold-value erreicht
  // double tstart=laser_sigma_t*sqrt(-2.0*log(tm2_laser_threshold));
  double tstart=tau_p * sqrt(-log(tm2_laser_threshold) / (4 * log(2)));

  laser_t_0 = tstart;
  laser_t_0 += 100e-15; // plus 100 fs puffer
#endif

#if defined(TMM) || defined(FDTD) || defined(LASER)
  error("TM2 does not work with TMM, FDTD or LASER. Only one of those options is allowed.");
#endif

  laser_spot_area=(fd_h.y*1e-10*fd_h.z*1e-10);

  if(myid==0){
    printf("*****************************************\n");
    printf("*        TRANSFER-MATRIX METHOD         *\n");
    printf("*****************************************\n");
    printf("lambda: %f nm\n", lambda*1.0e9);
    printf("theta: %f deg\n", tm2_theta);
    printf("tm2_absorption_threshold: %f\n", tm2_absorption_threshold);
    printf("tm2_laser_threshold: %f\n", tm2_laser_threshold);

    if(tm2_pol==1)
      printf("Polarization: S\n");
    if(tm2_pol==2)
      printf("Polarization: P\n");

    printf("I0: %.4e W/m^2\n", I0); // TODO: %f or %.4e
    printf("t0: %f fs\n", laser_t_0*1.0e15);
    printf("t1: %f fs\n", laser_t_1*1.0e15);
    printf("t_FWHM: %f fs, t1_FWHM: %f fs\n", tau_p*1.0e15, tau_p1*1.0e15);

#ifdef BALLISTIC
    printf("tau_relax: %f fs, lambda_ball: %f nm\n", tau_relax*1.0e15, lambda_ball*1.0e9);
#else
    if (lambda_ball>0)
      printf("tau_relax: %f fs, lambda_ball: %f nm\n", tau_relax*1.0e15, lambda_ball*1.0e9);
#endif

    printf("***************************\n");
  }

  return 0;
}


int do_tm2(real dt)
{
  int i,iglobal,j;

  I_t = exp(-4*log(2)*(tm2_time-laser_t_0)*(tm2_time-laser_t_0) / (tau_p*tau_p));

  if(steps<2)
  {
    return 0;
  }
 
  if(I_t<tm2_laser_threshold)
  {
    laser_active=false;
    return 0;
  }
  else
  {
    laser_active=true;
  }


  for(i=1;i<local_fd_dim.x-1;i++)
  {
#ifndef TTM1D
    iglobal = (i-1) + my_coord.x*(local_fd_dim.x-2); 
#else
    //bei TTM mit LB werden cpu-zuständigkeiten unabh. von position von links nach rechts
    //aufsteigend verteilt, sodass proc0 die zellen ganz links und proc "num_cpus-1" 
    //das rechte ende der simbox abarbeitet
    iglobal = (i-1) + myid*(local_fd_dim.x-2); 
#endif

    if(node.natoms < fd_min_atoms || node.dens < RHOMIN)
    {
      tm2_eps_real_arr_local[iglobal]=0.0;
      tm2_eps_imag_arr_local[iglobal]=0.0;
      tm2_active_cell_local[iglobal]=0;
      continue;
    }

#ifdef TFM
    tm2_get_epsilon(lambda, node.temp, node.md_temp,
                    node.Z, node.ne, &tm2_eps_real_arr_local[iglobal],
                    &tm2_eps_imag_arr_local[iglobal]);
#else
#ifdef TTT
    tm2_get_epsilon(lambda, node.temp, node.md_temp, &tm2_eps_real_arr_local[iglobal], &tm2_eps_imag_arr_local[iglobal], node.dens);
#endif // TTT
#endif // TFM
  }
  
  //MPI_Reduce(&tm2_alpha_arr_local[0], &tm2_alpha_arr_global[0],
  //           global_fd_dim.x, MPI_DOUBLE, MPI_SUM, 0, cpugrid);
  MPI_Allreduce(&tm2_alpha_arr_local[0], &tm2_alpha_arr_global[0], global_fd_dim.x, MPI_DOUBLE, MPI_SUM, cpugrid);
  MPI_Reduce(&tm2_active_cell_local[0], &tm2_active_cell_global[0],
             global_fd_dim.x, MPI_INT,MPI_SUM, 0, cpugrid);

  // Get reflectivity at room temperature and ambient density
  double eps_real_tmp, eps_imag_tmp;

#ifdef MATERIALAC
  float Z_tmp = MeanCharge(0.025851*11604.5, 2475, atomic_charge, atomic_weight[0], i, 1, 1);
  float ne_tmp = Z_tmp * 2475 / (atomic_weight[0] * AMU);
#else
  float Z_tmp = MeanCharge(0.025851*11604.5, 2475, atomic_charge, atomic_weight, i, 1, 1);
  float ne_tmp = Z_tmp * 2475 / (atomic_weight * AMU);
#endif // MATERIALAC

#ifdef TFM
  tm2_get_epsilon(lambda, 0.025851, 0.025851, Z_tmp, ne_tmp, &eps_real_tmp, &eps_imag_tmp);
#else
#ifdef TTT
#ifdef MATERIALCU
  tm2_get_epsilon(lambda, 0.025851, 0.025851, &eps_real_tmp, &eps_imag_tmp, 8628); // TODO: Density should be determined automatically
#else // MATERIALAL
  tm2_get_epsilon(lambda, 0.025851, 0.025851, &eps_real_tmp, &eps_imag_tmp, 2475); // TODO: Density should be determined automatically
#endif // MATERIALCU
#endif // TTT
#endif // TFM
  
  double n_tmp = sqrt((eps_real_tmp + sqrt(eps_real_tmp*eps_real_tmp + eps_imag_tmp*eps_imag_tmp)) / 2);
  double k_tmp = sqrt((-eps_real_tmp + sqrt(eps_real_tmp*eps_real_tmp + eps_imag_tmp*eps_imag_tmp)) / 2);
  double R_RT = ((n_tmp-1) * (n_tmp-1) + k_tmp*k_tmp) / ((n_tmp+1)*(n_tmp+1) + k_tmp*k_tmp);
  double alpha_tmp = 4 * M_PI * k_tmp / lambda;

  tm2_refl = R_RT; // TODO: Correct this

  // Find surface
  int last_empty_cell = -1; // if last_empty_cell==-1, all cells are filled, else last_empty_cell+1 is the surface cell
  for (int i=0; i<global_fd_dim.x-1; i++)
  {
    if (tm2_alpha_arr_global[i]==0.0) // Empty cell
    {
       last_empty_cell = i;
    }
    else
    {
      break; // There might be empty cells at the back of the sample which have to be ignored. Therefore, we end the search as soon as we've found some material
    }
  }

  for (int i=1; i<local_fd_dim.x-1; i++)
  {
#ifndef TTM1D
    iglobal = (i-1) + my_coord.x*(local_fd_dim.x-2); 
#else
    //bei TTM mit LB werden cpu-zuständigkeiten unabh. von position von links nach rechts
    //aufsteigend verteilt, sodass proc0 die zellen ganz links und proc "num_cpus-1" 
    //das rechte ende der simbox abarbeitet
    iglobal = (i-1) + myid*(local_fd_dim.x-2); 
#endif

    double alpha_sum = 0.0;
    for (int i2=0; i2<=iglobal; i2++) // Sum over all previous cells
    {
      alpha_sum += 1 / (1/(tm2_alpha_arr_global[i2]) + lambda_ball);
    }
    alpha_sum *= fd_h.x*1.0e-10; // Now it's the integral

    if (alpha_sum==0) // Empty cell
    {
      tm2_Qabs[iglobal] = 0.0;
    }
    else // Filled cell
    {
      tm2_Qabs[iglobal] = (1-R_RT) * 2 * sqrt(log(2)) * fluence * I_t * exp(-alpha_sum) / (sqrt(M_PI) * tau_p * (1/tm2_alpha_arr_global[iglobal] + lambda_ball));
    }
    tm2_Qabs[iglobal]*=6.3538562638e-26; // W/m^3 to imd-units

    node.source=node2.source=tm2_Qabs[iglobal];
  }

  return 0;
}


int tm2_get_epsilon_TFM(double lambda, double Te,double Ti, double Z, double Ne, real* re, real* img)
{
        Te*=11604.52500617; //ev->K
        Ti*=11604.52500617; //ev->K

        double complex un = 1;
        double complex im = I;

        double omega_las = (2.0 * M_PI * LIGHTSPEED / lambda);  // [1/s]
        double Ncr = ne_cr(omega_las);                  // [1/m3]
        double Ni = Ne / Z;                             // [1/m3]
        double EF = fermi_E(Ne);                        // [J]    
        double TF = 2.0 * EF / (3.0 * BOLTZMAN);        // [K] 
        double VF = sqrt(2 * EF / EMASS);               // [m/s] 
        double A1p=4.41;
        double A2p=0.8;
        double A3p=0.7;
        double A4p=0.2;



        int il, ir, id;

        for(id=0; lambda*1e6 > eps_bb_data[id][0]; id++){ // lambda [m]->[mu m]                 
                il = id;
                ir = id +1;
        }
        double nu_met=0;

        nu_met = numet(A1p, A2p, Te, Ti, TF);
        double nu_max=0;
        nu_max= numax(A3p, VF, Te, Ni);
        //if(isnan(nu_max)!=0)
        //  error("nu_max is NaN (getEps)");
        double nu_pl=0;
        nu_pl = nupl(omega_las, Z, Ni, Ne, Te);
        //if(isnan(nu_pl)!=0)
        //  error("nu_pl is NaN (getEps)");

        double ksi=0;
        ksi = 0.75 * sqrt(M_PI) * nu_pl / omega_las;


        double complex epsilon_bb = epsl_bb(il) + (epsl_bb(ir) - epsl_bb(il))/ (eps_bb_data[ir][0] - eps_bb_data[il][0])* (lambda*1e6 - eps_bb_data[il][0]);
        
#ifdef MATERIALAC
        int phase=((int) EOS_phase_from_r_ti(Ni*atomic_weight[0]*AMU,Ti));
#else
        int phase=((int) EOS_phase_from_r_ti(Ni*atomic_weight*AMU,Ti));
#endif // MATERIALAC

        if( fabs(phase)==3 || fabs(phase)==4 || fabs(phase)==5)
          epsilon_bb=0*I*0;


        double complex epsilon_met = epsilon_bb + un - Ne / Ncr / (un + im * MIN(nu_met, nu_max) / omega_las);
        double complex epsilon_pl  = un - Ne / Ncr * (tm2_K1(ksi) - im * nu_pl / omega_las * tm2_K2(ksi));
        double complex epsilon_wr = epsilon_pl + (epsilon_met - epsilon_pl) * exp(-A4p * Te / TF);

       *re=creal(epsilon_wr);
       *img=cimag(epsilon_wr);

//TESTCASE      
//Drude-Lorentz static
/*       
real epsinf=2.73;
real omega2= 2.2955e+16*2.2955e+16; //omega_plasma^2
real w2=omega_laser*omega_laser;
real gammap=1.1174e+15;
real omegap_L=7.6595e+15;
real omega0_L=2.4024e+15;
real gamma_L=4.5199e+14;
*/

//Drude-Lorentz dynamic, interpol
/*       
real epsinf=l1[i][j][1].DL[0];
real omega2=l1[i][j][1].DL[5]*l1[i][j][1].DL[5];
real w2=omega_laser*omega_laser;
real gammap=l1[i][j][1].DL[1];
real omegap_L=l1[i][j][1].DL[2];
real omega0_L=l1[i][j][1].DL[3];
real gamma_L=l1[i][j][1].DL[4];

double complex eps_dl= epsinf-omega2/(w2+I*gammap*omega_laser) + omegap_L*omegap_L/(omega0_L*omega0_L-w2-I*gamma_L*omega_laser);

*re= creal(eps_dl);
*img=cimag(eps_dl);   
*/
        return 1;
}

double complex epsl_bb(int bi){
        return eps_bb_data[bi][1] + I* eps_bb_data[bi][2];
}


#ifdef TTT
#ifdef MATERIALAL
int tm2_get_epsilon_Al(double lambda, double Te, double Ti, real* re, real* img, real rho)
{
  double E_Fermi = 11.1;  // in eV

  double eta = 0.75;
  double a0 = 0.49188e15;   // in 1/s
  double a1 = 1.2;
  double a2 = 0.8;
  double p1 = 0.75;
  double p2 = 2.25;

  double omega_rt_sol = 1.167022e14; // in 1/s
  double eta_sol = 0.5175;

  double omega_si_liq = 1.09084e15; // in 1/s
  double eta_liq = 0.4025;

  double omega_D = 16.4868e15;  // in 1/s

  double A_1 = 2.17020801;
  double A_2 = 3.95505075;
  double Omega_1 = 2.25028e15;  // in 1/s
  double Omega_2 = 2.26699e15;  // in 1/s
  double Gamma_1 = 0.268726e15; // in 1/s
  double Gamma_2 = 1.16111e15;  // in 1/s
  double phi_1 = -0.52279535;
  double phi_2 = 0.12358051;

  double x_rt = 0.991638025594; // 2.71213 / 2.735
  double x_m = 0.895795246801; // 2.45 / 2.735;

  double T_rt = 0.0258520402; // Room temperature in eV
  double T_m = 0.080439465748; // 933.5 / 11605 // Careful, this technically depends on the potential

  double x = 1; // TODO: rho / 2475;  // rho / rho_0 with rho_0 = 2475 kg/m^3
  double t = 6 * Te / (pow(x, 2/3) * E_Fermi);

  double omega_se = eta * pow(x, 2/3) * a0 * pow(t, 2) * (1 + a1 * pow(t, p1)) / (1 + a2 * pow(t, p2));

  // Determine whether material is solid or not via phase boundary lines
  int solid_bool = 1; // 1: Solid, 0: liquid // TODO: Proper evaluation of the state

  double omega_si;
  if (solid_bool == 1)
  {
    double y_bar_x = calc_y_bar(x);
    double y_bar_x_rt = calc_y_bar(x_rt);

    omega_si = omega_rt_sol * (1 + eta_sol * pow(x_rt, 2/3) * y_bar_x_rt * Ti / (pow(x, 2/3) * y_bar_x * T_rt));
  }
  else  // liquid
  {
    double beta = calc_beta(Ti);
    double w_Ti = calc_w(Ti);
    double w_Tm = calc_w(T_m);

    omega_si = eta_liq * omega_si_liq * pow(x, 2/3) * pow(x_m, beta) * w_Ti / (pow(x_m, 2/3) * pow(x, beta) * w_Tm);
  }

  double omega_eff = omega_se + omega_si;
  double omega = 2 * M_PI * LIGHTSPEED / lambda;  // in 1/s

  double complex im = I;

  double complex epsilon_Drude = 1 - x * omega_D * omega_D / (omega * omega + I * omega * omega_eff);
  double complex epsilon_DCP_1 = x * A_1 * Omega_1 * (exp(I * phi_1) / (Omega_1 - omega - I*Gamma_1)
                                                    + exp(-I * phi_1) / (Omega_1 + omega + I * Gamma_1));
  double complex epsilon_DCP_2 = x * A_2 * Omega_2 * (exp(I * phi_2) / (Omega_2 - omega - I*Gamma_2)
                                                    + exp(-I * phi_2) / (Omega_2 + omega + I * Gamma_2));

  double complex epsilon_tot = epsilon_Drude; // + epsilon_DCP_1 + epsilon_DCP_2;       // TODO: Add DCP

  *re = creal(epsilon_tot);
  *img = cimag(epsilon_tot);

  return 1;
}
#endif // MATERIALAL

#ifdef MATERIALCU
int tm2_get_epsilon_Cu(double lambda, double Te, double Ti, real* re, real* img, real rho)
{
  double E_Fermi = 9.2; // in eV

  // Parmeters from from Ren/Winter et al. for temperature- and density-independent optical properties
  double epsilon_infinity = 3.686;
  double omega_D = 1.34e16; // in rad/s, units don't matter here too much since we divide by the conversion factors
  double gamma = 0.99e14;
  double B_1 = 0.562;
  double B_2 = 27.36;
  double B_3 = 0.242;
  double phi_1 = -8.185;
  double phi_2 = 0.226;
  double phi_3 = -0.516;
  double Omega_1 = 3.205e15; // in rad/s
  double Omega_2 = 3.43e15;
  double Omega_3 = 7.33e15;
  double Gamma_1 = 0.404e15; // in rad/s
  double Gamma_2 = 0.77e16;
  double Gamma_3 = 1.12e15;

  // Extensions of optical properties by Winter et al.
  double T_rt = 0.0252488;  // 293 K in eV
  double T_m = 0.1006514714; // 1168.01 / 11604.5 // NOTE: This depends on the potential, this is the value of the unexcited ANA potential

  double C_dc = 1.79;
  double C_b = 1.031;
  double eta = 0.32;
  double beta_1 = 2.5e-2; // in 1/K
  double beta_2 = 3.5e-5; // in 1/K
  double n_at = 8.5e28; // in 1/m^3
  double m_opt = 2.63 * EMASS; // in kg

  double x = rho / 8628;  // rho / rho_rt with rho_rt = 8628 kg/m^3
  double t = 6 * Te / (x * E_Fermi);

  // Temperature and density-dependent Drude
  double omega_ss = 0.623 * x * t * t * (1 + 4.034 * pow(t, 0.7439)) * 1.0e15 / (1 + 1.2178 * pow(t, 2.0939)); // Also convert from 1/fs to 1/s
  double omega_sd = 5.32 * x * exp(-1.738 / t) * (1 + 1.64 * pow(t, -0.761)) * 1.0e15 / (1 + 0.0384 * pow(t, 0.673)); // Also convert from 1/fs to 1/s
  double omega_si = (C_dc + C_b * (1 - pow(MIN(Ti/T_m, 1), 0.5))) * Ti / 6.582119569e-16; // Last value is the reduced Planck's constant in eV * s
  double omega_eff = eta * (omega_ss + omega_sd) + omega_si;

  // Temperature-dependent collisions frequencies for DCP
  Gamma_1 *= (1 + beta_1 * (Ti - T_rt) * 11604.5);
  Gamma_2 *= (1 + beta_2 * (Ti - T_rt) * 11604.5);

  // Temperature and density-dependent plasma frequency
  // MYMOD: Simon, 17.12.2025
  double N_e;
  if (Te<=0.5)
  {
    N_e = 1.9;
  }
  else
  {
    N_e = 1.9 + 0.1485*(Te - 0.5);
  }
  // ENDOF MYMOD

  double omega_p = sqrt(N_e * n_at * ECHARGE * ECHARGE * x / (8.854188e-12 * m_opt));
  omega_D = omega_p;
  double omega = 2 * M_PI * LIGHTSPEED / lambda;  // in 1/s
  
  // MYMOD: Simon, 12.02.2026
  double epsilon_Drude_real = epsilon_infinity - pow(omega_D, 2) / (pow(omega, 2) + pow(omega_eff, 2));
  double epsilon_Drude_imaginary = pow(omega_D, 2) * omega_eff / (pow(omega, 3) + pow(omega_eff, 2) * omega);
  
  double epsilon_DCP_11_real = B_1 * Omega_1 * (cos(phi_1) * (Omega_1 - omega) - sin(phi_1) * Gamma_1) / (pow(Omega_1 - omega, 2) + pow(Gamma_1, 2));
  double epsilon_DCP_12_real = B_1 * Omega_1 * (cos(phi_1) * (Omega_1 + omega) - sin(phi_1) * Gamma_1) / (pow(Omega_1 + omega, 2) + pow(Gamma_1, 2));
  double epsilon_DCP_21_real = B_2 * Omega_2 * (cos(phi_2) * (Omega_2 - omega) - sin(phi_2) * Gamma_2) / (pow(Omega_2 - omega, 2) + pow(Gamma_2, 2));
  double epsilon_DCP_22_real = B_2 * Omega_2 * (cos(phi_2) * (Omega_2 + omega) - sin(phi_2) * Gamma_2) / (pow(Omega_2 + omega, 2) + pow(Gamma_2, 2));
  double epsilon_DCP_31_real = B_3 * Omega_3 * (cos(phi_3) * (Omega_3 - omega) - sin(phi_3) * Gamma_3) / (pow(Omega_3 - omega, 2) + pow(Gamma_3, 2));
  double epsilon_DCP_32_real = B_3 * Omega_3 * (cos(phi_3) * (Omega_3 + omega) - sin(phi_3) * Gamma_3) / (pow(Omega_3 + omega, 2) + pow(Gamma_3, 2));
  
  double epsilon_DCP_11_imaginary = B_1 * Omega_1 * (sin(phi_1) * (Omega_1 - omega) + cos(phi_1) * Gamma_1) / (pow(Omega_1 - omega, 2) + pow(Gamma_1, 2));
  double epsilon_DCP_12_imaginary = B_1 * Omega_1 * (-(sin(phi_1) * (Omega_1 + omega) + cos(phi_1) * Gamma_1)) / (pow(Omega_1 + omega, 2) + pow(Gamma_1, 2));
  double epsilon_DCP_21_imaginary = B_2 * Omega_2 * (sin(phi_2) * (Omega_2 - omega) + cos(phi_2) * Gamma_2) / (pow(Omega_2 - omega, 2) + pow(Gamma_2, 2));
  double epsilon_DCP_22_imaginary = B_2 * Omega_2 * (-(sin(phi_2) * (Omega_2 + omega) + cos(phi_2) * Gamma_2)) / (pow(Omega_2 + omega, 2) + pow(Gamma_2, 2));
  double epsilon_DCP_31_imaginary = B_3 * Omega_3 * (sin(phi_3) * (Omega_3 - omega) + cos(phi_3) * Gamma_3) / (pow(Omega_3 - omega, 2) + pow(Gamma_3, 2));
  double epsilon_DCP_32_imaginary = B_3 * Omega_3 * (-(sin(phi_3) * (Omega_3 + omega) + cos(phi_3) * Gamma_3)) / (pow(Omega_3 + omega, 2) + pow(Gamma_3, 2));

  double complex epsilon_tot_real = epsilon_Drude_real + x * (epsilon_DCP_11_real + epsilon_DCP_12_real + epsilon_DCP_21_real + \
                                                              epsilon_DCP_22_real + epsilon_DCP_31_real + epsilon_DCP_32_real);
  double complex epsilon_tot_imaginary = epsilon_Drude_imaginary + x * (epsilon_DCP_11_imaginary + epsilon_DCP_12_imaginary + epsilon_DCP_21_imaginary + \
                                                                        epsilon_DCP_22_imaginary + epsilon_DCP_31_imaginary + epsilon_DCP_32_imaginary);
  // ENDOF MYMOD

  *re = epsilon_tot_real;
  *img = epsilon_tot_imaginary;

  return 1;
}
#endif // MATERIALCU
#endif // TTT


#ifdef TFM
int tm2_read_arr(char* infilename,double ***arr,int cols,int *rowsout)
{

  FILE* infile;
  int i,j,ch;

  //first count lines
  if(myid==0)
  {
    infile = fopen(infilename, "r");
    if(infile==NULL)
      error_str("Cannot open file:%s\n", infilename);
    do
    {
       ch = fgetc(infile);
       if(ch == '\n')
         ++(*rowsout);
    } while (ch != EOF);
    rewind(infile);
  }
  MPI_Bcast(rowsout,1,MPI_INT,0,MPI_COMM_WORLD);

  //reallocate memory
  double *buf;     //1D-buffer for 2D-array
  //MPI_Alloc_mem((*rowsout)*2*sizeof(double), MPI_INFO_NULL, &buf);
  //MPI_Alloc_mem((*rowsout)*sizeof(**arr), MPI_INFO_NULL,arr);
  *arr = malloc(sizeof(**arr) * (*rowsout));
  buf= malloc(sizeof(double)*cols*(*rowsout));

  for (i = 0; i < (*rowsout); i++)
    //MPI_Alloc_mem(2*sizeof(***arr),MPI_INFO_NULL,(*arr)[i]);
    (*arr)[i] = malloc(sizeof ***arr * cols);
  // Now read
  if(myid==0)
  {
    for(i=0;i<(*rowsout);i++)
    {
      for(j=0;j<cols;j++)
      {
        fscanf(infile,"%lf",&(*arr)[i][j]);
        buf[i*cols+j]=(*arr)[i][j];
      }
    }
    fclose(infile);
  }
  MPI_Bcast(buf,(*rowsout)*cols,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //now reconstruct arr from buf
  for(i=0;i<(*rowsout);i++)
  {
    for(j=0;j<cols;j++)
    {
      (*arr)[i][j]=buf[i*cols+j];
    }
  }
  //MPI_Free_mem(buf);
  //for (i = 0; i < *rowsout; ++i)
  //    free((*arr)[i]);
  //free(*arr);
  free(buf);
  return 0;
}


/**************************************
*
* Gain semi empirical parameters
* K1 and K2 from table K12 and
* generate indicee il and ir
*
***************************************/
double tm2_K1(double nu) {
        double tmp;
        int i=0;

        if(nu < 1.0 / 30.0)
        return 1;
        else if (nu < 50) {

                i = (int) nu * 30 - 1;
                if(i<0) i =0;
                if(i > 1498) i = 1498;
                tmp=K12[i][0] + (30 * nu - i - 1) * (K12[i+1][0] - K12[i][0]);
        } else{
                tmp=315.0 / 8.0 / (nu * nu)  - 675675.0 / 64.0 / (nu * nu * nu * nu);
        }

        return tmp;
}
double tm2_K2(double nu) {
        double tmp;
        int i=0;

        if(nu < 1.0 / 30.0)
        return 1;
        else if (nu < 50) {
                i = (int) nu * 30 - 1;
                if(i<0) i=0;
                if(i > 1498) i = 1498;
                tmp=K12[i][1] + (30 * nu - i - 1) * (K12[i+1][1] - K12[i][1]);
        } else{
                tmp= 6.0 / (nu * nu)  - 720.0 / (nu * nu * nu * nu);
        }
        return tmp;
}
#endif // TFM
