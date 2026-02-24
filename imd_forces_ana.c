
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2024 Institute for Functional Matter and Quantum Technologies,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
*  All functions needed for ANA
*
******************************************************************************/

/******************************************************************************
*  $Revision$
*  $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/******************************************************************************
*
*  special version of do_forces for asymmetric core potentials
*  Needed for proper electron pressure
*
******************************************************************************/

#ifdef ASYMPOT

void do_forces(cell *p, cell *q, vektor pbc, real *Epot, real *Virial, 
               real *Vir_xx, real *Vir_yy, real *Vir_zz,
               real *Vir_yz, real *Vir_zx, real *Vir_xy)
{  
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real r2, rho_h;
  real tmp_virial;

  real pot_zwi, pot_grad;
  int inc = ntypes * ntypes;
  int jstart;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  
#ifdef ANAT
  int i_global, j_global;
  double t_el_ana_p, t_el_ana_q;
  double t_eff_p, t_eff_q;

  // Determine electron temperature in cell p
  if (p->n > 0)
  {
    i_global = (int) (ORT(p, 0, X)/fd_h.x); // Use particle 0 of second cell for determination of interaction/electron temperature
    // NOTE: If p->n == 0, this doesn't matter since no interaction is calculated and actio=reaction isn't used anyways
  }

  i_global = MIN(i_global, global_fd_dim.x-1);
  i_global = MAX(i_global, 0);
  t_eff_p = MIN(t_el_ana_arr_global[i_global], 1.2);    // ANA is only savely usable until ~1.2 eV

  // Determine electron temperature in cell q
  if (q->n > 0)
  {
    j_global = (int) (ORT(q, 0, X)/fd_h.x); // Use particle 0 of second cell for determination of interaction/electron temperature
    // NOTE: If q->n == 0, this doesn't matter since no interaction is calculated and actio=reaction isn't used anyways
  }

  j_global = MIN(j_global, global_fd_dim.x-1);
  j_global = MAX(j_global, 0);
  t_eff_q = MIN(t_el_ana_arr_global[j_global], 1.2);    // ANA is only savely usable until ~1.2 eV

  // If the electron temperature is the same in cell p and q, we only have to calculate the parameters once
  if (t_eff_p==t_eff_q){
    CALC_TE_PARAMS(t_eff_p);
  }
#endif // ANAT
  
  tmp_virial     = 0.0;
    
  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {
    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);

    /* for each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

      /* calculate distance */
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
      d.z = *qptr - tmp_d.z; ++qptr;
      r2    = SPROD(d,d);

      /* compute pair interactions, first on particle i */
      if (r2 <= ana_cutoff*ana_cutoff) {
#ifdef LOP
        // NOTE: Repeating the calculation of LOP isn't needed for particle j
        // since LOP is based on the structure, not on the electron temperature
        if (r2<=pow((sqrt(2)*a0/2 + a0)/2, 2)){
          // Nearest-neighbour vectors
          double nn_vectors[12][3] = {{a0*0.5, a0*0.5, 0.0},
                                      {-a0*0.5, a0*0.5, 0.0},
                                      {a0*0.5, -a0*0.5, 0.0},
                                      {-a0*0.5, -a0*0.5, 0.0},
                                      {0.0, a0*0.5, a0*0.5},
                                      {0.0, -a0*0.5, a0*0.5},
                                      {0.0, a0*0.5, -a0*0.5},
                                      {0.0, -a0*0.5, -a0*0.5},
                                      {a0*0.5, 0.0, a0*0.5},
                                      {-a0*0.5, 0.0, a0*0.5},
                                      {a0*0.5, 0.0, -a0*0.5},
                                      {-a0*0.5, 0.0, -a0*0.5}};

          double phi_tmp;
          double phi_min;
          
          for (int nn_vector_counter=0; nn_vector_counter<12; nn_vector_counter++){
            phi_tmp = pow(d.x - nn_vectors[nn_vector_counter][0], 2)
                    + pow(d.y - nn_vectors[nn_vector_counter][1], 2)
                    + pow(d.z - nn_vectors[nn_vector_counter][2], 2);

            // This is used to find the best agreement
            if (nn_vector_counter==0){
              phi_min = phi_tmp;
            }
            else{
              if (phi_tmp<phi_min){
                phi_min = phi_tmp;
              }
            }
          }

          // Assign best fitting value to arrays
          LOPI(p, i) += phi_min;
          LOPI(q, j) += phi_min;
          LOPN(p, i) += 1;
          LOPN(q, j) += 1;
        }
#endif // LOP

#ifdef TTM
        NUMNEIGHS(p,i)++;
        NUMNEIGHS(q,j)++;
#endif

#ifdef ANAT
        if (t_eff_p!=t_eff_q){
          CALC_TE_PARAMS(t_eff_q);
        }
#endif

        /* compute pair interactions */
        CALC_PAIR_ENERGY(pot_zwi, pot_grad);

        /* store force in temporary variable */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
        force.z = d.z * pot_grad;

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        *pfptr     += force.x; 
        *(++pfptr) += force.y; 
        *(++pfptr) += force.z; 

        /* the first half of the pot. energy of this bond */
        pot_zwi     *= 0.5;
        *Epot       += pot_zwi;
        POTENG(p,i) += pot_zwi;

        /* for the virial, we take the mean forces on the two particles */
        force.x       *= 0.5;
        force.y       *= 0.5;
        force.z       *= 0.5;
        pot_grad      *= 0.5;

        tmp_virial     -= r2 * pot_grad;

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) -= d.x * force.x;
          PRESSTENS(p,i,yy) -= d.y * force.y;
          PRESSTENS(p,i,zz) -= d.z * force.z;
          PRESSTENS(p,i,yz) -= d.y * force.z;
          PRESSTENS(p,i,zx) -= d.z * force.x;
          PRESSTENS(p,i,xy) -= d.x * force.y;
        }
#endif
      }

      /* compute pair interactions, now on particle j */
      if (r2 <= ana_cutoff*ana_cutoff) {
#ifdef ANAT
        if (t_eff_p!=t_eff_q){
          CALC_TE_PARAMS(t_eff_p);
        }
#endif
        /* compute pair interactions */
        CALC_PAIR_ENERGY(pot_zwi, pot_grad);

        /* store force in temporary variable */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
        force.z = d.z * pot_grad;

        /* accumulate forces */
        qfptr = &KRAFT(q,j,X);
        *qfptr     -= force.x; 
        *(++qfptr) -= force.y; 
        *(++qfptr) -= force.z; 

        /* the second half of the pot. energy of this bond */
        pot_zwi     *= 0.5;
        *Epot       += pot_zwi;
        POTENG(q,j) += pot_zwi;

        /* for the virial, we take the mean forces on the two particles */
        force.x       *= 0.5;
        force.y       *= 0.5;
        force.z       *= 0.5;
        pot_grad      *= 0.5;

        tmp_virial     -= r2 * pot_grad;

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(q,j,xx) -= d.x * force.x;
          PRESSTENS(q,j,yy) -= d.y * force.y;
          PRESSTENS(q,j,zz) -= d.z * force.z;
          PRESSTENS(q,j,yz) -= d.y * force.z;
          PRESSTENS(q,j,zx) -= d.z * force.x;
          PRESSTENS(q,j,xy) -= d.x * force.y;
        }
#endif
      }

      /* compute host electron density */
      if (r2 < ana_cutoff*ana_cutoff)  {
#ifdef ANAT
        if (t_eff_p!=t_eff_q){
          CALC_TE_PARAMS(t_eff_q);
        }
#endif
        CALC_RHO(rho_h);
        ANA_RHO(p,i) += rho_h; 

#ifdef ANAT
        if (t_eff_p!=t_eff_q){
          CALC_TE_PARAMS(t_eff_p);
        }
#endif
        CALC_RHO(rho_h);
        ANA_RHO(q,j) += rho_h; 
      }
    } /* for j */
  } /* for i */

  *Virial += tmp_virial;
}

#endif // ASYMPOT

/******************************************************************************
*
*  compute embedding energy and its derivative for all atoms
*
******************************************************************************/

void do_embedding_energy_ana(void)
{
  int k;
  
#ifdef ANAT
  int i_global;
  double t_el_ana, t_eff;
#endif  

  for (k=0; k<NCELLS; k++) {
    int  i, idummy=0;
    real pot;
    cell *p;
    p = CELLPTR(k);
    
#ifdef ANAT
    i_global=(int) (ORT(p,0,X)/fd_h.x);
    i_global=MIN(i_global,global_fd_dim.x-1);
    i_global=MAX(i_global,0);          
    t_el_ana=t_el_ana_arr_global[i_global];
    t_eff = MIN(t_el_ana, 1.2);    // ANA is only savely usable until ~1.2 eV

    CALC_TE_PARAMS(t_eff);
#endif

    for (i=0; i<p->n; i++) {
      if (ANA_RHO(p, i)==0.0){
        pot = 0.0;
      }
      else{
        CALC_EMBEDDING_ENERGY(pot, ANA_DF(p,i), ANA_RHO(p,i));
      
        POTENG(p,i)    += pot;
        tot_pot_energy += pot;
      }
    }
  }
}

/******************************************************************************
*
*  second force loop, calculates the force and the energy 
*  caused by the embedding electron density
*  uses Phi(r2), Rho(r2), F(rho) and its derivatives
*
******************************************************************************/

void do_forces_ana(cell *p, cell *q, vektor pbc, real *Virial, 
                    real *Vir_xx, real *Vir_yy, real *Vir_zz,
                    real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,k,same_cell;
  vektor d, tmp_d, force;
  real r2;
  int  idummy=0;
  int  jstart;
  int  inc=ntypes*ntypes;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  real tmp_virial=0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif
  real ana_force, rho_i_strich, rho_j_strich;
 
#ifdef ANAT
  int i_global, j_global;
  double t_el_ana_p, t_el_ana_q;
  double t_eff_p, t_eff_q;

  // Determine electron temperature in cell p
  if (p->n > 0)
  {
    i_global = (int) (ORT(p, 0, X)/fd_h.x); // Use particle 0 of second cell for determination of interaction/electron temperature
    // NOTE: If p->n == 0, this doesn't matter since no interaction is calculated and actio=reaction isn't used anyways
  }

  i_global = MIN(i_global, global_fd_dim.x-1);
  i_global = MAX(i_global, 0);
  t_eff_p = MIN(t_el_ana_arr_global[i_global], 1.2);    // ANA is only savely usable until ~1.2 eV

  // Determine electron temperature in cell q
  if (q->n > 0)
  {
    j_global = (int) (ORT(q, 0, X)/fd_h.x); // Use particle 0 of second cell for determination of interaction/electron temperature
    // NOTE: If q->n == 0, this doesn't matter since no interaction is calculated and actio=reaction isn't used anyways
  }

  j_global = MIN(j_global, global_fd_dim.x-1);
  j_global = MAX(j_global, 0);
  t_eff_q = MIN(t_el_ana_arr_global[j_global], 1.2);    // ANA is only savely usable until ~1.2 eV

  // If the electron temperature is the same in cell p and q, we only have to calculate the parameters once
  if (t_eff_p==t_eff_q){
    CALC_TE_PARAMS(t_eff_p);
  }
#endif // ANAT

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {
    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    same_cell = ((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0));

    jstart = (same_cell ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);

    /* for each atom in neighbouring cell */
    for (j=jstart; j<q->n; ++j) {

      /* calculate distance */ 
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
      d.z = *qptr - tmp_d.z; ++qptr;
      r2    = SPROD(d,d);
      
      if (r2 < ana_cutoff * ana_cutoff) {
        /* take care: particle i gets its rho from particle j. */

#ifdef ANAT
        if (t_eff_p!=t_eff_q){
          CALC_TE_PARAMS(t_eff_q);
        }
#endif

        /* rho_strich_i(r_ij) */
        CALC_RHO_DERIV(rho_i_strich);

        /* rho_strich_j(r_ij) */
#if defined(ANAT)
        if (t_eff_p==t_eff_q) {
          rho_j_strich = rho_i_strich;
        } else {
          CALC_TE_PARAMS(t_eff_p);
          CALC_RHO_DERIV(rho_j_strich);
        }
#else
        rho_j_strich = rho_i_strich;
#endif
        
        /* put together (dF_i and dF_j are by 0.5 too big) */
        ana_force = 0.5 * (ANA_DF(p,i)*rho_j_strich + ANA_DF(q,j)*rho_i_strich);

        /* store force in temporary variable */
        force.x = d.x * ana_force;
        force.y = d.y * ana_force;
        force.z = d.z * ana_force;

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        qfptr = &KRAFT(q,j,X);
        *pfptr     += force.x; 
        *qfptr     -= force.x; 
        *(++pfptr) += force.y; 
        *(++qfptr) -= force.y; 
        *(++pfptr) += force.z; 
        *(++qfptr) -= force.z; 

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
        tmp_vir_vect.z -= d.z * force.z;
#else
        tmp_virial     -= r2  * ana_force;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          /* avoid double counting of the virial */
          force.x *= 0.5;
          force.y *= 0.5;
          force.z *= 0.5;
 
          PRESSTENS(p,i,xx) -= d.x * force.x;
          PRESSTENS(p,i,yy) -= d.y * force.y;
          PRESSTENS(p,i,zz) -= d.z * force.z;
          PRESSTENS(p,i,yz) -= d.y * force.z;
          PRESSTENS(p,i,zx) -= d.z * force.x;
          PRESSTENS(p,i,xy) -= d.x * force.y;

          PRESSTENS(q,j,xx) -= d.x * force.x;
          PRESSTENS(q,j,yy) -= d.y * force.y;
          PRESSTENS(q,j,zz) -= d.z * force.z;
          PRESSTENS(q,j,yz) -= d.y * force.z;
          PRESSTENS(q,j,zx) -= d.z * force.x;
          PRESSTENS(q,j,xy) -= d.x * force.y;
	}
#endif

      } /* if in the cutoff range */

    } /* for j */
  } /* for i */

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif 

} /* do_forces_ana */

/******************************************************************************
*
*  Correct the potential energies to a fixed reference energy
*
******************************************************************************/

#ifdef ANAT

void correct_energy_anat(void)
{
  int k;
  int i_global;
  double t_el_ana, t_eff;

  for (k=0; k<NCELLS; k++) {
    int  i, idummy=0;
    real pot;
    cell *p;
    p = CELLPTR(k);

    i_global=(int) (ORT(p,0,X)/fd_h.x);
    i_global=MIN(i_global,global_fd_dim.x-1);
    i_global=MAX(i_global,0);
    t_el_ana=t_el_ana_arr_global[i_global];
    t_eff = MIN(t_el_ana, 1.2);    // ANA is only savely usable until ~1.2 eV

    // Single-atom energies
    double E_single_RT = -0.23359741 - 0.36458457*0.025851 - 2.36130430*pow(0.025851, 2);
    double E_single_T = -0.23359741 - 0.36458457*t_eff - 2.36130430*pow(t_eff, 2);

    for (i=0; i<p->n; i++) {
      // Correct potential energy to fixed single-atom energy (RT)
      POTENG(p, i) += E_single_T - E_single_RT;
    }
  }
}

#endif

