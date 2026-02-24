
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
*  Macros for potential table access. We use macros in order to enforce
*  inlining, which is hard to achieve with functions. We need inlining
*  because these access functions are called very frequently.
*
*  $Revision$
*  $Date$
*
******************************************************************************/

/* 4-point (cubic), spline, or 3-point (quadratic) interpolation of 
   function tables */
#if   defined(FOURPOINT)
#define   PAIR_INT   PAIR_INT3
#define   VAL_FUNC   VAL_FUNC3
#define DERIV_FUNC DERIV_FUNC3
#elif defined(SPLINE)
#define   PAIR_INT   PAIR_INT_SP
#define   VAL_FUNC   VAL_FUNC_SP
#define DERIV_FUNC DERIV_FUNC_SP
#else
#define   PAIR_INT   PAIR_INT2
#define   VAL_FUNC   VAL_FUNC2
#define DERIV_FUNC DERIV_FUNC2
#endif

/* compensate for non-standard cast in icc (which is faster) */
/* this works only for a positive argument */
#ifdef RCD
#define POS_TRUNC(x) ((int)(-0.5+(x)))
#else
#define POS_TRUNC(x) ((int)(x))
#endif

/*****************************************************************************
*
*  Evaluate harmonic crystal
*
******************************************************************************/

#define ATOM_INT_EC(pot, grad, p_typ, r2)		        	\
{									\
  pot += 0.5 * spring_rate[p_typ] * r2;			        	\
  grad += spring_rate[p_typ];				        	\
}

/*****************************************************************************
*
*  Evaluate Lennard-Jones potential
*
******************************************************************************/

#define PAIR_INT_LJ(pot, grad, p_typ, q_typ, r2)			\
{									\
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;				\
									\
  sig_d_rad2  = lj_sigma[p_typ][q_typ] * lj_sigma[p_typ][q_typ] / (r2);	\
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;			\
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;				\
									\
  pot += lj_epsilon[p_typ][q_typ] * ( sig_d_rad12 - 2.0 * sig_d_rad6 )	\
    - lj_shift[p_typ][q_typ];						\
  grad += - 12.0 * lj_epsilon[p_typ][q_typ] / (r2)			\
    * ( sig_d_rad12 - sig_d_rad6 );					\
}

#ifdef VEC

#define PAIR_INT_LJ_VEC(pot, grad, col, r2)		                   \
{									   \
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;				   \
									   \
  sig_d_rad2  = lj_sigma2_vec[col] / (r2);                                 \
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;			   \
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;				   \
									   \
  pot  = lj_epsilon_vec[col] * (sig_d_rad12 - 2.0 * sig_d_rad6);           \
  grad = - 12.0 * lj_epsilon_vec[col] * (sig_d_rad12 - sig_d_rad6) / (r2); \
}

/*****************************************************************************
*
*  Evaluate Lennard-Jones potential
*
******************************************************************************/

#define PAIR_INT_LJG(pot, grad, p_typ, q_typ, r2)			\
{									\
  real sig_d_rad2, sig_d_rad6, sig_d_rad12, expo, dr;			\
									\
  sig_d_rad2  = lj_sigma[p_typ][q_typ] * lj_sigma[p_typ][q_typ] / (r2);	\
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;			\
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;				\
                                                                        \
  dr = ( sqrt(r2)-ljg_r0[p_typ][q_typ] ) / ljg_sig[p_typ][q_typ];       \
  expo = exp ( - 0.5 * dr * dr );					\
									\
  pot += lj_epsilon[p_typ][q_typ] * ( sig_d_rad12 - 2.0 * sig_d_rad6 )	\
    - lj_shift[p_typ][q_typ]                                            \
    - ljg_eps[p_typ][q_typ] * expo;					\
  grad += - 12.0 * lj_epsilon[p_typ][q_typ] / (r2)			\
    * ( sig_d_rad12 - sig_d_rad6 )                                      \
    - ljg_eps[p_typ][q_typ] * dr * expo / ljg_sig[p_typ][q_typ];	\

/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT_VEC(pot, grd, pt, col, inc, r2)                             \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  long kk;                                                                   \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  kk    = (long) (r2a);                                                      \
  chi   = r2a - kk;                                                          \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, kk, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
  grd = 2 * istep * (dv + (chi - 0.5) * d2v);                                \
}

#define VAL_FUNC_VEC(pot, pt, col, inc, r2)                                  \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  long kk;                                                                   \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  kk    = (long) (r2a);                                                      \
  chi   = r2a - kk;                                                          \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, kk, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential value */                                                      \
  pot = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
}

#define DERIV_FUNC_VEC(grd, pt, col, inc, r2)                                \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  long kk;                                                                   \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  kk    = (long) (r2a);                                                      \
  chi   = r2a - kk;                                                          \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, kk, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* twice the derivative */                                                 \
  grd = 2 * istep * (dv + (chi - 0.5) * d2v);                                \
}

#endif

/*****************************************************************************
*
*  Evaluate Morse potential
*
******************************************************************************/

#define PAIR_INT_MORSE(pot, grad, p_typ, q_typ, r2)			    \
{									    \
  real r, exppot, cexppot;						    \
									    \
  r       = sqrt((r2));							    \
  exppot  = exp( - morse_alpha[p_typ][q_typ]				    \
		 * ( r - morse_sigma[p_typ][q_typ] ) );			    \
  cexppot = 1.0 - exppot;						    \
									    \
  pot  += morse_epsilon[p_typ][q_typ] * ( cexppot * cexppot - 1.0 )	    \
    - morse_shift[p_typ][q_typ];					    \
  grad += 2.0 * morse_alpha[p_typ][q_typ] * morse_epsilon[p_typ][q_typ] / r \
    * exppot * cexppot;							    \
}

/*****************************************************************************
*
*  Evaluate Buckingham potential 
*
******************************************************************************/

#define PAIR_INT_BUCK(pot, grad, p_typ, q_typ, r2)			   \
{									   \
  real rinv, rinv2, powpot, exppot, invs2;				   \
									   \
  rinv   = buck_sigma[p_typ][q_typ] / sqrt((r2));			   \
  rinv2  = rinv * rinv;							   \
  powpot  = buck_c[p_typ][q_typ] * rinv2 * rinv2 * rinv2;		   \
  exppot = buck_a[p_typ][q_typ] * exp ( - 1.0 / rinv );			   \
  invs2   = 1.0 / ( buck_sigma[p_typ][q_typ] * buck_sigma[p_typ][q_typ] ); \
									   \
  pot += exppot - powpot - buck_shift[p_typ][q_typ];			   \
  grad += ( - exppot * rinv + 6 * powpot * rinv2 ) * invs2;		   \
}

/*****************************************************************************
*
*  Evaluate pair potential for Keating potential 
*
******************************************************************************/

#define PAIR_INT_KEATING(pot, grad, p_typ, q_typ, r2)                         \
{                                                                             \
  real tmp_pot, tmp;                                                          \
                                                                              \
  tmp_pot = 3.0 * keat_alpha[p_typ][q_typ]                                    \
          / ( 8.0 * keat_d[p_typ][q_typ] * keat_d[p_typ][q_typ] );            \
  tmp = r2 - keat_d[p_typ][q_typ] * keat_d[p_typ][q_typ];                     \
                                                                              \
  pot = tmp_pot * tmp * tmp;                                                  \
  grad = 4.0 * tmp_pot * tmp;                                                 \
}

/*****************************************************************************
*
*  Evaluate pair potential for Stillinger-Weber potential 
*
******************************************************************************/

#define PAIR_INT_STIWEB(pot, grad, p_typ, q_typ, r2)                          \
{                                                                             \
  real radius, phi_r, phi_a, inv_c, inv_r, f_cut;                             \
                                                                              \
  radius = sqrt((r2));                                                        \
  phi_r  = sw_a[p_typ][q_typ] * pow( radius, - sw_p[p_typ][q_typ] );          \
  phi_a  = - sw_b[p_typ][q_typ] * pow( radius, - sw_q[p_typ][q_typ] );        \
  inv_c  = 1.0 / ( (radius) - sw_a1[p_typ][q_typ] );                          \
  inv_r  = 1.0 / (radius);                                                    \
  f_cut  = exp( sw_de[p_typ][q_typ] * inv_c );                                \
                                                                              \
  pot  = ( phi_r + phi_a ) * f_cut;                                           \
  grad = ( - pot * sw_de[p_typ][q_typ] * inv_c * inv_c                        \
         - f_cut * inv_r * ( sw_p[p_typ][q_typ] * phi_r                       \
                           + sw_q[p_typ][q_typ] * phi_a ) ) * inv_r;          \
                                                                              \
}


/*****************************************************************************
*
*  Evaluate potential table with linear interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT_LIN(pot, grad, pt, col, inc, r2, is_short)                  \
{                                                                            \
  real r2a, chi, *t;                                                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  r2a   = r2a *  (pt).invstep[col];                                          \
  k     = POS_TRUNC(r2a);                                                    \
  chi   = r2a - k;                                                           \
                                                                             \
  /* potential and twice the derivative */                                   \
  t    = (pt).table[col];                                                    \
  pot  = t[2*k  ]*(1.0-chi) + t[2*k+2]*chi;                                  \
  grad = t[2*k+1]*(1.0-chi) + t[2*k+3]*chi;                                  \
}

/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT2(pot, grad, pt, col, inc, r2, is_short)                     \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                        \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/*****************************************************************************
*
*  Evaluate potential table with cubic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT3(pot, grad, pt, col, inc, r2, is_short)                     \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;                                \
  real fac0, fac1, fac2, fac3, dfac0, dfac1, dfac2, dfac3;                   \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = MAX( POS_TRUNC(r2a), 1 );                                          \
  /* k     = MIN( MAX( POS_TRUNC(r2a), 1 ), (pt).len[col]-3 ); */            \
  chi   = r2a - k;                                                           \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                           \
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);                             \
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);                           \
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);                                   \
                                                                             \
  /* factors for the interpolation of the derivative */                      \
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);                              \
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);                              \
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);                              \
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);                                    \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k-1, (col), (pt).maxsteps, (inc));                \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* potential energy */                                                     \
  pot = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;                       \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);    \
}

/*****************************************************************************
*
*  Evaluate potential table with spline interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT_SP(pot, grad, pt, col, inc, r2, is_short)                   \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  step  = (pt).step[col];                                                    \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-2 ); */                      \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k     = k * (inc) + (col);                                                 \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k    += (inc);                                                             \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = a * p1 + b * p2 + (a * a2 * d21 + b * b2 * d22) * st6 * step;       \
  grad = 2*((p2 - p1) * istep + ((3*b2 + 2) * d22 - (3*a2 + 2) * d21) * st6);\
}

/*****************************************************************************
*
*  Evaluate tabulated function with quadratic interpolation. 
*
******************************************************************************/

#define VAL_FUNC2(val, pt, col, inc, r2, is_short)                           \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* the function value */                                                   \
  val = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
}

/*****************************************************************************
*
*  Evaluate tabulated function with cubic interpolation. 
*
******************************************************************************/

#define VAL_FUNC3(val, pt, col, inc, r2, is_short)                           \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3;                                      \
  real fac0, fac1, fac2, fac3, *ptr;                                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = MAX( POS_TRUNC(r2a), 1 );                                          \
  /* k     = MIN( MAX( POS_TRUNC(r2a), 1 ), (pt).len[col]-3 ); */            \
  chi   = r2a - k;                                                           \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                           \
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);                             \
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);                           \
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);                                   \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k-1, (col), (pt).maxsteps, (inc));                \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* the function value */                                                   \
  val = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;                       \
}

/*****************************************************************************
*
*  Evaluate tabulated function with spline interpolation. 
*
******************************************************************************/

#define VAL_FUNC_SP(val, pt, col, inc, r2, is_short)                         \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  step  = (pt).step[col];                                                    \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-2 ); */                      \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k     = k * (inc) + (col);                                                 \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k    += (inc);                                                             \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
                                                                             \
  /* the function value */                                                   \
  val  = a * p1 + b * p2 + (a * a2 * d21 + b * b2 * d22) * st6 * step;       \
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with quadratic interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

#define DERIV_FUNC2(grad, pt, col, inc, r2, is_short)                        \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a   = 0;                                                               \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with cubic interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

#define DERIV_FUNC3(grad, pt, col, inc, r2, is_short)                        \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;                                \
  real dfac0, dfac1, dfac2, dfac3;                                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = MAX( POS_TRUNC(r2a), 1 );                                          \
  /* k     = MIN( MAX( POS_TRUNC(r2a), 1 ), (pt).len[col]-3 ); */            \
  chi   = r2a - k;                                                           \
                                                                             \
  /* factors for the interpolation of the 1. derivative */                   \
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);                              \
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);                              \
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);                              \
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);                                    \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k-1, (col), (pt).maxsteps, (inc));                \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);    \
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with spline interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

#define DERIV_FUNC_SP(grad, pt, col, inc, r2, is_short)                      \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  step  = (pt).step[col];                                                    \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-2 ); */                      \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k     = k * (inc) + (col);                                                 \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k    += (inc);                                                             \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2*((p2 - p1) * istep + ((3*b2 + 2) * d22 - (3*a2 + 2) * d21) * st6);\
}

#ifdef ANA

/******************************************************************************
*
*  compute pair energy and its derivative for all atoms
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*
******************************************************************************/

#define CALC_PAIR_ENERGY(pot_zwi, pot_grad)                                   \
{                                                                             \
  double r = sqrt(r2);                                                        \
  double x = ( r - ana_cutoff ) / sc_h;                                       \
  double x3 = x * x * x;                                                      \
  double x4 = x3 * x;                                                         \
  double exp_morse = exp(- morse_a * ( r - morse_r0 ));                       \
  double one_exp_morse = 1 - exp_morse;                                       \
  double one_plus_x4 = 1 + x4;                                                \
  double one_over_one_plus_x4 = 1 / one_plus_x4;                              \
  double x3_times_one_over_one_plus_x4 = x3 * one_over_one_plus_x4;           \
  double x4_times_one_over_one_plus_x4 = x3_times_one_over_one_plus_x4 * x;   \
                                                                              \
  pot_zwi = morse_D                                                           \
            * ( one_exp_morse * one_exp_morse - 1)                            \
            * x4_times_one_over_one_plus_x4;                                  \
  pot_grad = 2 * morse_D * morse_a * exp_morse * one_exp_morse                \
                         * x4_times_one_over_one_plus_x4 / r                  \
             + pot_zwi * ( 4 * x3_times_one_over_one_plus_x4                  \
                           - 4 * x3_times_one_over_one_plus_x4                \
                               * x4_times_one_over_one_plus_x4 )              \
                       / ( r * sc_h );                                        \
}

/******************************************************************************
*
*  compute embedding energy and its derivative for one pair of atoms
*
******************************************************************************/

#define CALC_EMBEDDING_ENERGY(embed_e_c, embed_e_c_grad, rho_c)               \
{                                                                             \
  double log_rho = log(rho_c);                                                \
  double rho_to_gamma_1 = pow(rho_c, bjs_gamma - 1);                          \
                                                                              \
  embed_e_c = bjs_f0 * ( 1 - bjs_gamma * log_rho ) * rho_to_gamma_1 * rho_c   \
              + bjs_f1 * rho_c;                                               \
  embed_e_c_grad = 2 * ( - bjs_f0 * bjs_gamma * bjs_gamma * log_rho           \
                                  * rho_to_gamma_1 + bjs_f1 );                \
}

#define CALC_EMBEDDING_ENERGY_(embed_e_c, embed_e_c_grad, rho_c)              \
{                                                                             \
  embed_e_c = 0.0;                                                            \
  embed_e_c_grad = 0.0;                                                       \
}

/******************************************************************************
*
*  compute host electron density
*
******************************************************************************/

#define CALC_RHO(rho_c)                                                       \
{                                                                             \
  double r = sqrt(r2);                                                        \
  double x = ( r - ana_cutoff ) / sc_h;                                       \
  double x4 = x * x * x * x;                                                  \
                                                                              \
  rho_c = ( 1 + csw_a1 * cos(csw_alpha * r) + csw_a2 * sin(csw_alpha * r) )   \
          * ( x4 / ( 1 + x4) ) / pow(r, csw_beta);                            \
}

#define CALC_RHO_(rho_c)                                                      \
{                                                                             \
  rho_c = 0.0;                                                                \
}

/******************************************************************************
*
*  compute host electron density derivative
*
******************************************************************************/

#define CALC_RHO_DERIV(drho_c)                                                \
{                                                                             \
  double r = sqrt(r2);                                                        \
  double x = ( r - ana_cutoff ) / sc_h;                                       \
  double x3 = x * x * x;                                                      \
  double x4 = x3 * x;                                                         \
  double cos_ar = cos(csw_alpha * r);                                         \
  double sin_ar = sin(csw_alpha * r);                                         \
  double r_to_beta = pow(r, csw_beta);                                        \
  double r_to_beta_plus_one = r_to_beta * r;                                  \
  double one_plus_x4 = 1 + x4;                                                \
  double one_over_one_plus_x4 = 1 / one_plus_x4;                              \
  double x3_times_one_over_one_plus_x4 = x3 * one_over_one_plus_x4;           \
  double x4_times_one_over_one_plus_x4 = x * x3_times_one_over_one_plus_x4;   \
                                                                              \
  double rho_c_no_sc = ( 1 + csw_a1 * cos_ar + csw_a2 * sin_ar ) / r_to_beta; \
  drho_c = ( ( csw_a2 * csw_alpha * cos_ar - csw_a1 * csw_alpha * sin_ar )    \
             / r_to_beta_plus_one                                             \
             - csw_beta * rho_c_no_sc / r2 ) * x4_times_one_over_one_plus_x4  \
           + rho_c_no_sc * ( 4 * x3_times_one_over_one_plus_x4                \
                             - 4 * x3_times_one_over_one_plus_x4              \
                                 * x4_times_one_over_one_plus_x4 )            \
                           / ( r * sc_h );                                    \
}

#define CALC_RHO_DERIV_(drho_c)                                               \
{                                                                             \
  drho_c = 0.0;                                                               \
}

/******************************************************************************
*
*  compute electron temperature-dependent parameters for ANA
*
******************************************************************************/

#define CALC_TE_PARAMS(E_el)                                                  \
{                                                                             \
  double D_e_0 = 0.13176571148031033;                                         \
  double D_e_1 = 0.01754114533241918;                                         \
  double D_e_2 = -0.09007937629407885;                                        \
  double D_e_3 = 0.0205682667797867;                                          \
                                                                              \
  double a_0 = 1.9415027938258338;                                            \
  double a_1 = 0.008527493893889335;                                          \
  double a_2 = -0.04182841074885612;                                          \
  double a_3 = 0.054999099334988456;                                          \
                                                                              \
  double r_e_0 = 2.676898870235321;                                           \
  double r_e_1 = -0.01704117223687041;                                        \
  double r_e_2 = 0.09799670817457108;                                         \
  double r_e_3 = 0.029374520746831487;                                        \
                                                                              \
  double a1_0 = -0.414829562655917;                                           \
  double a1_1 = 0.41420238402813064;                                          \
  double a1_2 = -1.0767485471886642;                                          \
  double a1_3 = 1.1707157794466336;                                           \
                                                                              \
  double a2_0 = -0.4087441558313223;                                          \
  double a2_1 = -0.16050673266534451;                                         \
  double a2_2 = -0.3362426127432448;                                          \
  double a2_3 = 0.2760693880261744;                                           \
                                                                              \
  double alpha_0 = 2.0530885578802596;                                        \
  double alpha_1 = 0.11030379504738536;                                       \
  double alpha_2 = -0.1544388140691187;                                       \
  double alpha_3 = 0.30432216428156617;                                       \
                                                                              \
  double beta_0 = 3.3298091580028575;                                         \
  double beta_1 = 0.1932335043978177;                                         \
  double beta_2 = -0.3969886857297274;                                        \
  double beta_3 = 0.27890286039150314;                                        \
                                                                              \
  double F_0_0 = -2.538473239507734;                                          \
  double F_0_1 = 0.9037531239692363;                                          \
  double F_0_2 = 0.5050490366157505;                                          \
  double F_0_3 = 0.06982969253330407;                                         \
                                                                              \
  double gamma_0 = 0.7934921852755044;                                        \
  double gamma_1 = 0.41031019891622533;                                       \
  double gamma_2 = -0.8092964040574772;                                       \
  double gamma_3 = 0.6770372352871102;                                        \
                                                                              \
  double F_1_0 = -0.012107002017232884;                                       \
  double F_1_1 = -0.013463096415215283;                                       \
  double F_1_2 = 0.02690409718255482;                                         \
  double F_1_3 = -0.02088764700760113;                                        \
                                                                              \
  morse_D = D_e_0 + D_e_1 * E_el + D_e_2 * E_el * E_el + D_e_3 * E_el * E_el * E_el;              \
  morse_a = a_0 + a_1 * E_el + a_2 * E_el * E_el + a_3 * E_el * E_el * E_el;                      \
  morse_r0 = r_e_0 + r_e_1 * E_el + r_e_2 * E_el * E_el + r_e_3 * E_el * E_el * E_el;             \
  csw_a1 = a1_0 + a1_1 * E_el + a1_2 * E_el * E_el + a1_3 * E_el * E_el * E_el;                   \
  csw_a2 = a2_0 + a2_1 * E_el + a2_2 * E_el * E_el + a2_3 * E_el * E_el * E_el;                   \
  csw_alpha = alpha_0 + alpha_1 * E_el + alpha_2 * E_el * E_el + alpha_3 * E_el * E_el * E_el;    \
  csw_beta = beta_0 + beta_1 * E_el + beta_2 * E_el * E_el + beta_3 * E_el * E_el * E_el;         \
  bjs_f0 = F_0_0 + F_0_1 * E_el + F_0_2 * E_el * E_el + F_0_3 * E_el * E_el * E_el;               \
  bjs_gamma = gamma_0 + gamma_1 * E_el + gamma_2 * E_el * E_el + gamma_3 * E_el * E_el * E_el;    \
  bjs_f1 = F_1_0 + F_1_1 * E_el + F_1_2 * E_el * E_el + F_1_3 * E_el * E_el * E_el;               \
}

#endif
