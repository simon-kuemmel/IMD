#include "imd.h"
#include <complex.h>
/*
	Transfer-Matrix Methode zum Lösen der Helmholtz-Wellengleichung.
*/

// ACHTUNG: Falls Qabs=NaN --> Probe wahrscheinlich zu kurz!

 
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
#define tmm_get_epsilon tmm_get_epsilon_TFM
#endif
#ifdef TTT
#ifdef MATERIALAL // Al
#define tmm_get_epsilon tmm_get_epsilon_Al
#define calc_y_bar calc_y_bar_Al
#define calc_beta calc_beta_Al
#define calc_w calc_w_Al
#endif // MATERIALAL
#ifdef MATERIALCU // Cu
#define tmm_get_epsilon tmm_get_epsilon_Cu
#define calc_y_bar calc_y_bar_Cu
#define calc_beta calc_beta_Cu
#define calc_w calc_w_Cu
#endif // MATERIALCU
#endif

/*2x2 Matrix multiplication*/
void matmul(double complex *a,double complex *b,double complex *p)
{
        p[0]=a[0]*b[0]+a[2]*b[1];
        p[1]=a[1]*b[0]+a[3]*b[1];
        p[2]=a[0]*b[2]+a[2]*b[3];
        p[3]=a[1]*b[2]+a[3]*b[3];
}
/* 2x2 matrix and 2x1 vector multiplication */
void matmul2(double complex *mat,double complex *vec,double complex *out)
{
        out[0]=mat[0]*vec[0]+mat[2]*vec[1];
        out[1]=mat[1]*vec[0]+mat[3]*vec[1];
}
double Runge5(double delta, double complex kl, double complex epsl, double complex Bplus, double complex Bminus) 
{
  double dx = 1;
  double cur_xpos = 0;
  double result = 0;
  const double errval = 1e-5;
  double k1, k3, k4, k5;


  while(cur_xpos < 1 && dx > 1e-5) {

        k1 = dx / 3.0 * EE(cur_xpos, delta, kl, epsl, Bplus, Bminus);
        k3 = dx / 3.0 * EE(cur_xpos + dx / 3.0, delta, kl, epsl, Bplus, Bminus);
        k4 = dx / 3.0 * EE(cur_xpos + 0.5 * dx, delta, kl, epsl, Bplus, Bminus);
        k5 = dx / 3.0 * EE(cur_xpos + dx, delta, kl, epsl, Bplus, Bminus);
        double ERR = k1 - 4.5 * k3 + 4 * k4 - 0.5 * k5;

        if(ERR < 5.0 / 32 * errval) {
                result += 0.5 * (k1 + 4 * k4 + k5);
                cur_xpos += dx;
                dx *= 1.1;
        } else if (ERR  > 5 * errval) {
                dx *= 0.5;
                continue;
        } else {
                result += 0.5 * (k1 + 4 * k4 + k5);
                cur_xpos += dx;
        }

        if(cur_xpos + dx > 1) dx = 1 - cur_xpos;

        }
  return result;
}


double EE(double z, double delta, double complex kl, double complex epsl, double complex Bplus, double complex Bminus) {
        double complex un=1+I*0;
        double complex im=0+I*1;

        //const double sin_theta = sin(m_theta * PI / 180.0);
        double phi1   = creal(delta * z * kl * im);
        double phi2   = cimag(delta * z * kl * im);
        double complex eiphi=exp(phi1)*cos(phi2)+I*exp(phi1)*sin(phi2);

        double complex Ez, Ey;
        /*
        if (m_polarization == 'p') {
                Ez = sin_theta / epsl * (Bplus * eiphi + Bminus / eiphi );
                Ey = - cmplroot(epsl - sin_theta * sin_theta) / epsl * (Bplus * eiphi - Bminus / eiphi);
        } else {
                Ez = (Bplus * eiphi + Bminus / eiphi);
                Ey = 0;
        }
        */

        Ez = (Bplus * eiphi + Bminus / eiphi); //<--eigentlich ist das Ey.
        Ey = 0;				     

        return creal(Ez)*creal(Ez)+cimag(Ez)*cimag(Ez) + creal(Ey)*creal(Ey)+cimag(Ey)*cimag(Ey);
}

int tmm_init()
{ 
  //alloc eps arr
  // tmm_eps_real_arr_local = (real*) malloc(global_fd_dim.x*sizeof(real));
  // tmm_eps_imag_arr_local = (real*) malloc(global_fd_dim.x*sizeof(real));

  // tmm_eps_real_arr_global= (real*) malloc(global_fd_dim.x*sizeof(real));
  // tmm_eps_imag_arr_global= (real*) malloc(global_fd_dim.x*sizeof(real));

  alloc1darr(real,tmm_eps_real_arr_global,global_fd_dim.x);
  alloc1darr(real,tmm_eps_imag_arr_global,global_fd_dim.x);

  alloc1darr(real,tmm_eps_real_arr_local,global_fd_dim.x);
  alloc1darr(real,tmm_eps_imag_arr_local,global_fd_dim.x);  
  
#ifdef TFM
  tmm_read_arr(eps_filename,&eps_bb_data,3,&eps_bb_rows);
  tmm_read_arr(k12_filename,&K12,2,&K12_rows); //K1 und K1 num.eval.integrale
#endif

  // km=(double complex*) malloc(global_fd_dim.x*sizeof(double complex));

  alloc1darr(double complex,km,global_fd_dim.x);

  // tmm_active_cell_local=(int*) malloc(global_fd_dim.x*sizeof(int));
  // tmm_active_cell_global=(int*) malloc(global_fd_dim.x*sizeof(int));

  alloc1darr(int,tmm_active_cell_global,global_fd_dim.x);
  alloc1darr(int,tmm_active_cell_local,global_fd_dim.x);

  // tmm_Qabs =      (double*) calloc(global_fd_dim.x,sizeof(double));
  // tmm_Qabs_scat = (double*) calloc(global_fd_dim.x,sizeof(double));

  alloc1darr(double,tmm_Qabs,global_fd_dim.x);
  alloc1darr(double,tmm_Qabs_scat,global_fd_dim.x);
  int i;

  tmm_dt=timestep*10.18/1e15; // in sek.
  for(i=0;i<global_fd_dim.x;++i)
  {
    tmm_eps_real_arr_local[i]=tmm_eps_real_arr_global[i]=0.0;
    tmm_eps_imag_arr_local[i]=tmm_eps_imag_arr_global[i]=0.0;   
    km[i]=0.0+I*0.0;
    tmm_active_cell_local[i]=tmm_active_cell_global[i]=0;

    tmm_Qabs[i]=0.0;
    tmm_Qabs_scat[i]=0.0;
  }

  real freq=c0/lambda;
  omega_laser=2.0*M_PI*freq;
  k0=2.0*M_PI/lambda;
  k0/=1e10; //weil fd_h.x in Angstrom 


#ifdef TMM_t0_suggest
  // verschiebe t0 sodasss laser E-feld gerade threshold-value erreicht
  //double tstart=laser_sigma_t*sqrt(-2.0*log(tmm_laser_threshold));  
  double tstart=tau_p * sqrt(-log(tmm_laser_threshold) / (4 * log(2)));
  laser_t_0 = tstart;
  laser_t_0 += 100e-15; // plus 100 fs puffer
#endif


#if defined(FDTD) || defined(LASER)
    error("TMM does not work with FDTD or LASER. Only one of those options is allowed.");
#endif

    laser_spot_area=(fd_h.y*1e-10*fd_h.z*1e-10);
        if(myid==0){
          printf("*****************************************\n");
          printf("*        TRANSFER-MATRIX METHOD         *\n");
          printf("*****************************************\n");
          printf("lambda:%.4e\n",lambda);
          printf("theta:%f (deg)\n", tmm_theta);
          printf("tmm_absorption_threshold:%f\n",tmm_absorption_threshold);
          printf("tmm_laser_threshold:%f\n",tmm_laser_threshold);
          if(tmm_pol==1)
            printf("Polarization: S\n");
          if(tmm_pol==2)
	        printf("Polarization: P\n");
          printf("I0:%.4e W/m^2\n",I0);
          printf("t0:%.4e\n",laser_t_0);
          printf("t1:%.4e\n",laser_t_1);
          printf("t_FWHM:%.4e, t1_FWHM:%.4e\n", tau_p, tau_p1);          
          printf("***************************\n");
	}

  return 0;
}
// **********************************************************************************************************************
int do_tmm(real dt)
{  
  int i,iglobal,j;
  double Imp=3.769911184307751e+02; // vacuum impedanz
  I_t=I0*exp(-4*log(2)*pow(tmm_time-laser_t_0,2)/(tau_p*tau_p));
  I_t+=I0*exp(-4*log(2)*pow(tmm_time-laser_t_1,2)/(tau_p1*tau_p1));
  
  if(steps<2) 
  {
    return 0;
  }
  
  if (I_t<I0*tmm_laser_threshold)
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

   //if(l1[i][1][1].natoms<fd_min_atoms)
  if(node.natoms < fd_min_atoms || node.dens < RHOMIN)
   {
     tmm_eps_real_arr_local[iglobal]=0.0;
     tmm_eps_imag_arr_local[iglobal]=0.0;
     tmm_active_cell_local[iglobal]=0; 
     continue;
   } 
   //tmm_eps_real_arr_local[iglobal]=ttm_get_epsilon();
   //tmm_eps_imag_arr_local[iglobal]=i


   // tmm_get_epsilon(lambda,0.0258, 0.0258,
   //                2.5, 2e29, &tmm_eps_real_arr_local[iglobal],
   //                &tmm_eps_imag_arr_local[iglobal]);

#ifdef TFM
   tmm_get_epsilon(lambda, node.temp, node.md_temp,
                  node.Z, node.ne, &tmm_eps_real_arr_local[iglobal],
                  &tmm_eps_imag_arr_local[iglobal]);
#endif
#ifdef TTT
    tmm_get_epsilon(lambda, node.temp, node.md_temp, &tmm_eps_real_arr_local[iglobal], &tmm_eps_imag_arr_local[iglobal], node.dens);
#endif

    tmm_active_cell_local[iglobal]=1; 
 } 

 MPI_Reduce(&tmm_eps_real_arr_local[0], &tmm_eps_real_arr_global[0],
  	    global_fd_dim.x, MPI_DOUBLE,MPI_SUM, 0, cpugrid);
 MPI_Reduce(&tmm_eps_imag_arr_local[0], &tmm_eps_imag_arr_global[0],
            global_fd_dim.x, MPI_DOUBLE,MPI_SUM, 0, cpugrid);
 MPI_Reduce(&tmm_active_cell_local[0], &tmm_active_cell_global[0],
            global_fd_dim.x, MPI_INT,MPI_SUM, 0, cpugrid);

 //KORREKTUR-LOOP
 for(i=1;i<local_fd_dim.x-1;i++)
 {
#ifndef TTM1D
   iglobal = (i-1) + my_coord.x*(local_fd_dim.x-2);
#else
   iglobal = (i-1) + myid*(local_fd_dim.x-2);
#endif
   if (tmm_eps_real_arr_global[iglobal]==0) tmm_eps_real_arr_global[iglobal]=1.0; //VACUUM hat epsilon_real=1
 }
 // 
 //  
 // 
  double complex Cinv[4];
  double complex Pmat[4];
  double complex tmpmat[4];
  double complex tmpmat2[4];
  double complex Emat[4];

  double complex kl,kr;
  double phi1,phi2;
  double complex eiphi;

  double complex BR,BT,B0;
  int ecut=0;
  real ecut_thresh=exp(-tmm_absorption_threshold);

  BR=BT=0.0+I*0;
  B0=1.0+I*0;

  real dx=fd_h.x;
//////////////////////
//1st loop: calc km
//////////////////////
if(myid==0)
{
  km[0]=k0;
  for(i=1;i<global_fd_dim.x-1;i++)
  {
   if(tmm_active_cell_global[i]==0)
	   km[i]=k0;
   else
	   km[i]=k0*csqrt(tmm_eps_real_arr_global[i]+I*tmm_eps_imag_arr_global[i]);
 }
}
///////////////////////
//2nd loop: find ecut
/////////////////////
if(myid==0)
{
  Emat[0]=1.0+I*0;
  Emat[1]=0;
  Emat[2]=0;
  Emat[3]=1.0+I*0.0;


  int tooshort=1;
  for(i=0;i<global_fd_dim.x-1;i++)
  { 
    //km[i]=k0*csqrt(tmm_eps_real_arr_global[i]+I*tmm_eps_imag_arr_global[i]);
    if(i==0) dx=1e10;else dx=fd_h.x;
    kl=km[i];
    kr=km[i+1];
    phi1=creal(dx*kl*I);
    phi2=cimag(dx*kl*I);
    eiphi=exp(phi1)*cos(phi2)+I*exp(phi1)*sin(phi2);

    //Pmat
    Pmat[0]=eiphi;
    Pmat[2]=1.0/eiphi;
    Pmat[1]=eiphi;
    Pmat[3]=-Pmat[2];

    //Cinv
    Cinv[0]=0.5;
    Cinv[2]=0.5*kl/kr;
    Cinv[1]=0.5;
    Cinv[3]=-Cinv[2];  

    //E=(Cinv*P) * E
    matmul(Cinv,Pmat,tmpmat); 
    matmul(tmpmat,Emat,tmpmat2);
    for(j=0;j<4;j++)
	     Emat[j]=tmpmat2[j];
    
    BR=-Emat[1]*B0 / Emat[3];
    BT=Emat[0]*B0+ Emat[2]*BR;

    if(creal(BT)*creal(BT)+cimag(BT)*cimag(BT)< ecut_thresh)
    {
        tooshort=0;
	      ecut=i+1;
	      break;
    }
  }

  if(tooshort==1)
  {
    int l;
    for(l=0;l<global_fd_dim.x-1;l++)
    {
      printf("km:%f+i*%f\n", creal(km[l]),cimag(km[l]));
    }
    error("Sample is too short for TMM.Consider reducing tmm_threshold.");
    MPI_Abort(cpugrid,0);
  }

  //tans und refl
  double complex tmpr0=BR/B0;
  tmm_refl=creal(tmpr0)*creal(tmpr0)+cimag(tmpr0)*cimag(tmpr0);
  tmm_trans=creal(BT)*creal(BT)+cimag(BT)*cimag(BT);
  tmm_abs=1.0-tmm_refl-tmm_trans;
}
////////////////
// 3rd loop for Fvec
/////////////////
if(myid==0)
{
  double complex* Bplus= (double complex*) calloc(ecut+1,sizeof(double complex));
  double complex* Bminus=(double complex*) calloc(ecut+1,sizeof(double complex));
  double complex Fvec[2];
  double complex tmpvec[2];

  Fvec[0]=B0;
  Fvec[1]=BR;
  Bplus[0]=B0;
  Bminus[0]=BR;

  for(i=0;i<ecut;i++)
  {
    if(i==0) dx=1e10;else dx=fd_h.x;
    kl=km[i];
    kr=km[i+1];
    phi1=creal(dx*kl*I);
    phi2=cimag(dx*kl*I);
    eiphi=exp(phi1)*cos(phi2)+I*exp(phi1)*sin(phi2);

    //Pmat
    Pmat[0]=eiphi;
    Pmat[2]=1.0/eiphi;
    Pmat[1]=eiphi;
    Pmat[3]=-Pmat[2];

    //Cinv
    Cinv[0]=0.5;
    Cinv[2]=0.5*kl/kr;
    Cinv[1]=0.5;
    Cinv[3]=-Cinv[2];


    //F=(C*P)*F;
    matmul(Cinv,Pmat,tmpmat);
    matmul(tmpmat,Fvec,tmpvec);
    Fvec[0]=tmpvec[0]; Fvec[1]=tmpvec[1];

    Bplus[i+1]=Fvec[0];
    Bminus[i+1]=Fvec[1];
  }
//////////////////////
// 4th loop for Qabs
//////////////////////  
  for(i=0;i<ecut;i++)
  {
    if(i==0) dx=1e10;else dx=fd_h.x;

    kl=km[i];
    double complex eps=tmm_eps_real_arr_global[i]+I*tmm_eps_imag_arr_global[i];
    tmm_Qabs[i]=I_t*k0*tmm_eps_imag_arr_global[i]*Runge5(dx,kl,eps,Bplus[i],Bminus[i]);

/*
    phi1   = creal(dx * kl * I);
    phi2   = cimag(dx * kl * I);
    eiphi=exp(phi1)*cos(phi2)+I*exp(phi1)*sin(phi2);
    Ez = (Bplus[i] * eiphi + Bminus[i] / eiphi); // Ez ist complex!
    tmm_Qabs[i]=I_t*k0*tmm_eps_imag_arr_global[i]*(creal(Ez)*creal(Ez)+cimag(Ez)*cimag(Ez));
*/
    tmm_Qabs[i]*=1e10; // W/m^2/Angstrom to W/m^3
    tmm_Qabs[i]*=6.3538562638e-26; // W/m^3 to imd-units   
  }
  free(Bplus);
  free(Bminus);
} // if myid==0
///////////////////////
//  NOW COMM.
///////////////////////
MPI_Scatter(&tmm_Qabs[0], local_fd_dim.x-2, MPI_DOUBLE, &tmm_Qabs_scat[0], local_fd_dim.x-2, MPI_DOUBLE, 0, cpugrid);
for(i=1;i<local_fd_dim.x-1;i++)
{
  node.source=node2.source=tmm_Qabs_scat[i-1];
}
 return 0;
}

/*
int tmm_get_epsilon(real* re, real* img)
{ 
  real omegap=2.160294047429742e+16;
  real gammap=1e14; //Drude

gammap=1.7927354e14; //Drude-Lorentz
omegap=2.0312604e16;

  real gamma2=gammap*gammap;
  real omega2=omegap*omegap;
  real w2=omega_laser*omega_laser;
  real w4=w2*w2;
  //double complex eps=1-(omegap*omegap)/(1+omega_laser*omega_laser+I*omega_laser*gammap);
  real tmpre,tmpimg;

  //Drude only
  //tmpre=1.0-(omega2*w2)/(gamma2*w2+w4);
  //tmpimg=gammap*omega2/(gamma2*omega_laser+w2*omega_laser);


  //Drude-Lorentz
  real omegap_L=9.6929256835e15;
  real gamma_L=9.84485242e14;
  real omega0_L=2.47640578e15;
  real epsinf=2.73;
  double complex eps_dl= epsinf-omega2/(w2+I*gammap*omega_laser) + omegap_L*omegap_L/(omega0_L*omega0_L-w2-I*gamma_L*omega_laser);
//eps_dl=-61.312323928473738 +41.843997569604973*I;
eps_dl=-64.199494418202008 +81.919352477839695*I; //TESTCASE, sharp, Gamma_L halbiert, lambda=761 nm (==Omega0_L)
eps_dl=-73.655835558841446 +81.545418969195765*I; // testcase, sharp, omegap_L*2, omegap_D/2, lambda=650 nm

  tmpimg=cimag(eps_dl);
  tmpre=creal(eps_dl);

  //tmpre*=eps0;
  //tmpimg*=eps0;

  *re=tmpre;
  *img=tmpimg;
  //eps*=eps0;
  return 0; //alles ok
}
*/

int tmm_get_epsilon_TFM(double lambda, double Te,double Ti, double Z, double Ne, real* re, real* img)
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
        int phase=((int) EOS_phase_from_r_ti(Ni*atomic_weight*AMU,Ti));

        if( fabs(phase)==3 || fabs(phase)==4 || fabs(phase)==5)
          epsilon_bb=0*I*0;


        double complex epsilon_met = epsilon_bb + un - Ne / Ncr / (un + im * MIN(nu_met, nu_max) / omega_las);
        double complex epsilon_pl  = un - Ne / Ncr * (tmm_K1(ksi) - im * nu_pl / omega_las * tmm_K2(ksi));
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
int tmm_get_epsilon_Al(double lambda, double Te, double Ti, real* re, real* img, real rho)
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

  double x = 1;	// TODO: rho / 2475;  // rho / rho_0 with rho_0 = 2475 kg/m^3
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

  double complex epsilon_tot = epsilon_Drude; // + epsilon_DCP_1 + epsilon_DCP_2;	// TODO: Add DCP

  *re = creal(epsilon_tot);
  *img = cimag(epsilon_tot);

  return 1;
}
#endif // MATERIALAL

#ifdef MATERIALCU
int tmm_get_epsilon_Cu(double lambda, double Te, double Ti, real* re, real* img, real rho)
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


int tmm_read_arr(char* infilename,double ***arr,int cols,int *rowsout)
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
double tmm_K1(double nu) {
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
double tmm_K2(double nu) {
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

