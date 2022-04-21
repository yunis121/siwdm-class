/** @file background.c Documented background module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondence between redshift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantities {B} can be found as a function of conformal
 * time by integrating the background equations.
 *
 * 3. some other quantities that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be inferred analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantities {A} as a function of quantities {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantities {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most precise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other quantity, we can call background_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background quantities, we should call first
 * background_tau_of_z() and then background_at_tau().
 *
 *
 * In order to save time, background_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantities, or rarely useful
 * quantities). Each line in the interpolation table is a vector whose
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantities like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantities needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning
 * -# background_at_tau(), background_tau_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

/* EDIT: Added Self Interacting Dark Matter Module - Rafael Yunis - 2020
 *
 * The only changes made here to the overall functioning of the module is to add the 
 * calculation (interpolation from table) of the relaxation time.
 * 
 * In addition, all functions regarding inizialitation of the SIWDM module reside in this file,
 * which include things like redefining the integration rules, determining what type of 
 * cosmological scenario takes place, etc are invluded.
 *
 * Inizialization of the SIWDM modlue is better understood together with the file input.c, 
 * which includes some of the inizialization calls to these functions
 */


#include "background.h"

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_tau(
                      struct background *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < pba->tau_table[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table[0]);

  class_test(tau > pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dtau2_table,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table,
                                                        pba->bt_size,
                                                        pba->background_table,
                                                        pba->d2background_dtau2_table,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;
}

/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int background_tau_of_z(
                        struct background *pba,
                        double z,
                        double * tau
                        ) {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Background quantities at given \f$ a \f$.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the beginning of this file). In extended
 * cosmological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param pvecback_B    Input: vector containing all {B} type quantities (scale factor, ...)
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assumed to be already allocated)
 * @return the error status
 */

int background_functions(
                         struct background *pba,
                         double * pvecback_B, /* Vector containing all {B} quantities. */
                         short return_format,
                         double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                         ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
  /* scale factor relative to scale factor today */
  double a_rel;
  /* background ncdm quantities */
  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  /* index for n_ncdm species */
  int n_ncdm;
  /* fluid's time-dependent equation of state parameter */
  double w_fld, dw_over_da, integral_fld;
  /* scale factor */
  double a;
  /* scalar field quantities */
  double phi, phi_prime;

  /* SI quantities */
  //number of si species
  int n_si_ncdm;
  // relaxation time
  double reltime;

  /** - initialize local variables */
  a = pvecback_B[pba->index_bi_a];
  rho_tot = 0.;
  p_tot = 0.;
  rho_r=0.;
  rho_m=0.;
  a_rel = a / pba->a_today;

  class_test(a_rel <= 0.,
             pba->error_message,
             "a = %e instead of strictly positive",a_rel);

  /** - pass value of \f$ a\f$ to output */
  pvecback[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }

  /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    /* Pass value of rho_dcdm to output */
    pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
    rho_tot += pvecback[pba->index_bg_rho_dcdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_dcdm];
  }

  /* dr */
  if (pba->has_dr == _TRUE_) {
    /* Pass value of rho_dr to output */
    pvecback[pba->index_bg_rho_dr] = pvecback_B[pba->index_bi_rho_dr];
    rho_tot += pvecback[pba->index_bg_rho_dr];
    p_tot += (1./3.)*pvecback[pba->index_bg_rho_dr];
    rho_r += pvecback[pba->index_bg_rho_dr];
  }

  /* Scalar field */
  if (pba->has_scf == _TRUE_) {
    phi = pvecback_B[pba->index_bi_phi_scf];
    phi_prime = pvecback_B[pba->index_bi_phi_prime_scf];
    pvecback[pba->index_bg_phi_scf] = phi; // value of the scalar field phi
    pvecback[pba->index_bg_phi_prime_scf] = phi_prime; // value of the scalar field phi derivative wrt conformal time
    pvecback[pba->index_bg_V_scf] = V_scf(pba,phi); //V_scf(pba,phi); //write here potential as function of phi
    pvecback[pba->index_bg_dV_scf] = dV_scf(pba,phi); // dV_scf(pba,phi); //potential' as function of phi
    pvecback[pba->index_bg_ddV_scf] = ddV_scf(pba,phi); // ddV_scf(pba,phi); //potential'' as function of phi
    pvecback[pba->index_bg_rho_scf] = (phi_prime*phi_prime/(2*a*a) + V_scf(pba,phi))/3.; // energy of the scalar field. The field units are set automatically by setting the initial conditions
    pvecback[pba->index_bg_p_scf] =(phi_prime*phi_prime/(2*a*a) - V_scf(pba,phi))/3.; // pressure of the scalar field
    rho_tot += pvecback[pba->index_bg_rho_scf];
    p_tot += pvecback[pba->index_bg_p_scf];
    //divide relativistic & nonrelativistic (not very meaningful for oscillatory models)
    rho_r += 3.*pvecback[pba->index_bg_p_scf]; //field pressure contributes radiation
    rho_m += pvecback[pba->index_bg_rho_scf] - 3.* pvecback[pba->index_bg_p_scf]; //the rest contributes matter
    //printf(" a= %e, Omega_scf = %f, \n ",a_rel, pvecback[pba->index_bg_rho_scf]/rho_tot );
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {

    /* Loop over species: */
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

      /* function returning background ncdm[n_ncdm] quantities (only
         those for which non-NULL pointers are passed) */
      class_call(background_ncdm_momenta(
                                         pba,
                                         n_ncdm,
                                         pba->q_ncdm_bg[n_ncdm],
                                         pba->w_ncdm_bg[n_ncdm],
                                         pba->q_size_ncdm_bg[n_ncdm],
                                         pba->M_ncdm[n_ncdm],
                                         pba->factor_ncdm[n_ncdm],
                                         1./a_rel-1.,
                                         NULL,
                                         &rho_ncdm,
                                         &p_ncdm,
                                         NULL,
                                         &pseudo_p_ncdm),
                 pba->error_message,
                 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;

      /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;

      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
         to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    rho_tot += pvecback[pba->index_bg_rho_lambda];
    p_tot -= pvecback[pba->index_bg_rho_lambda];
  }

  /* fluid with w(a) and constant cs2 */
  if (pba->has_fld == _TRUE_) {

    /* get rho_fld from vector of integrated variables */
    pvecback[pba->index_bg_rho_fld] = pvecback_B[pba->index_bi_rho_fld];

    /* get w_fld from dedicated function */
    class_call(background_w_fld(pba,a,&w_fld,&dw_over_da,&integral_fld), pba->error_message, pba->error_message);
    pvecback[pba->index_bg_w_fld] = w_fld;

    // Obsolete: at the beginning, we had here the analytic integral solution corresponding to the case w=w0+w1(1-a/a0):
    // pvecback[pba->index_bg_rho_fld] = pba->Omega0_fld * pow(pba->H0,2) / pow(a_rel,3.*(1.+pba->w0_fld+pba->wa_fld)) * exp(3.*pba->wa_fld*(a_rel-1.));
    // But now everthing is integrated numerically for a given w_fld(a) defined in the function background_w_fld.

    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += w_fld * pvecback[pba->index_bg_rho_fld];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_ur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_ur];
    rho_r += pvecback[pba->index_bg_rho_ur];
  }

  /** - compute expansion rate H from Friedmann equation: this is the
      only place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of \f$ [3c^2/8\pi G] \f$, ie
      \f$ \rho_{class} = [8 \pi G \rho_{physical} / 3 c^2]\f$ */
  pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);

  /** - compute derivative of H with respect to conformal time */
  pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_tot;

  /* SI */
  /* Calculate relaxation time using background_si_ncdm_reltime routine, store in pvecback*/

  for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

    int n_si_ncdm = 0;

    if (pba->ncdm_si_type[n_ncdm] >= 1){

      class_call(background_si_ncdm_reltime(pba, &reltime, 1./a_rel-1., n_ncdm)
        ,pba->error_message
        ,pba->error_message);

      n_si_ncdm = pba->ncdm_si_index[n_ncdm];

      pvecback[pba->index_bg_taurel_si_ncdm1 + n_si_ncdm ] = reltime;

      //printf("taurel_index = %d, reltime=%e \n", pba->index_bg_taurel_si_ncdm1, reltime );

    }
    
  }

  /** - compute other quantities in the exhaustive, redundant format */
  if (return_format == pba->long_info) {

    /** - compute critical density */
    pvecback[pba->index_bg_rho_crit] = rho_tot-pba->K/a/a;
    class_test(pvecback[pba->index_bg_rho_crit] <= 0.,
               pba->error_message,
               "rho_crit = %e instead of strictly positive",pvecback[pba->index_bg_rho_crit]);

    /** - compute Omega_m */
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_tot;

    /* one can put other variables here */
    /*  */
    /*  */

    // Additional info output regarding si_ncdm
    if (pba->has_si_ncdm == _TRUE_){
      if (pba->background_verbose>5){
        for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){
          n_si_ncdm = pba->ncdm_si_index[n_ncdm];
          printf("background_functions called with long_info. Printing taurel info: T = %e, z = %e, taurel = %e\n",
            ( pba->T_cmb * pba->T_ncdm[n_ncdm] * (1./a_rel) ),
            1./a_rel-1.,
            pvecback[pba->index_bg_taurel_si_ncdm1+n_si_ncdm] );
        }
      }
    }

  }

  return _SUCCESS_;

}

/**
 * Single place where the fluid equation of state is
 * defined. Parameters of the function are passed through the
 * background structure. Generalisation to arbitrary functions should
 * be simple.
 *
 * @param pba            Input: pointer to background structure
 * @param a              Input: current value of scale factor
 * @param w_fld          Output: equation of state parameter w_fld(a)
 * @param dw_over_da_fld Output: function dw_fld/da
 * @param integral_fld   Output: function \f$ \int_{a}^{a_0} da 3(1+w_{fld})/a \f$
 * @return the error status
 */

int background_w_fld(
                     struct background * pba,
                     double a,
                     double * w_fld,
                     double * dw_over_da_fld,
                     double * integral_fld) {

  double Omega_ede = 0.;
  double dOmega_ede_over_da = 0.;
  double d2Omega_ede_over_da2 = 0.;
  double a_eq, Omega_r, Omega_m;

  /** - first, define the function w(a) */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *w_fld = pba->w0_fld + pba->wa_fld * (1. - a / pba->a_today);
    break;
  case EDE:
    // Omega_ede(a) taken from eq. (10) in 1706.00730
    Omega_ede = (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))
      /(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      + pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld));

    // d Omega_ede / d a taken analytically from the above
    dOmega_ede_over_da = - pba->Omega_EDE* 3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.)/(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      - (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))*(1.-pba->Omega0_fld)*3.*pba->w0_fld*pow(a,3.*pba->w0_fld-1.)/pow(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld),2)
      + pba->Omega_EDE*3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.);

    // find a_equality (needed because EDE tracks first radiation, then matter)
    Omega_r = pba->Omega0_g * (1. + 3.046 * 7./8.*pow(4./11.,4./3.)); // assumes LambdaCDM + eventually massive neutrinos so light that they are relativistic at equality; needs to be generalised later on.
    Omega_m = pba->Omega0_b;
    if (pba->has_cdm == _TRUE_) Omega_m += pba->Omega0_cdm;
    if (pba->has_dcdm == _TRUE_)
        class_stop(pba->error_message,"Early Dark Energy not compatible with decaying Dark Matter because we omitted to code the calculation of a_eq in that case, but it would not be difficult to add it if necessary, should be a matter of 5 minutes");
    a_eq = Omega_r/Omega_m; // assumes a flat universe with a=1 today
    class_stop(pba->error_message,"a_eq = %e, z_eq =%e\n",a_eq,1./a_eq-1.);

    // w_ede(a) taken from eq. (11) in 1706.00730
    *w_fld = - dOmega_ede_over_da*a/Omega_ede/3./(1.-Omega_ede)+a_eq/3./(a+a_eq);
    break;
  }


  /** - then, give the corresponding analytic derivative dw/da (used
        by perturbation equations; we could compute it numerically,
        but with a loss of precision; as long as there is a simple
        analytic expression of the derivative of the previous
        function, let's use it! */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *dw_over_da_fld = - pba->wa_fld / pba->a_today;
    break;
  case EDE:
    d2Omega_ede_over_da2 = 0.;
    *dw_over_da_fld = - d2Omega_ede_over_da2*a/3./(1.-Omega_ede)/Omega_ede
      - dOmega_ede_over_da/3./(1.-Omega_ede)/Omega_ede
      + dOmega_ede_over_da*dOmega_ede_over_da*a/3./(1.-Omega_ede)/(1.-Omega_ede)/Omega_ede
      + a_eq/3./(a+a_eq)/(a+a_eq);
    break;
  }

  /** - finally, give the analytic solution of the following integral:
        \f$ \int_{a}^{a0} da 3(1+w_{fld})/a \f$. This is used in only
        one place, in the initial conditions for the background, and
        with a=a_ini. If your w(a) does not lead to a simple analytic
        solution of this integral, no worry: instead of writing
        something here, the best would then be to leave it equal to
        zero, and then in background_initial_conditions() you should
        implement a numerical calculation of this integral only for
        a=a_ini, using for instance Romberg integration. It should be
        fast, simple, and accurate enough. */
  *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(pba->a_today/a) + pba->wa_fld*(a/pba->a_today-1.));

  /** note: of course you can generalise these formulas to anything,
      defining new parameters pba->w..._fld. Just remember that so
      far, HyRec explicitely assumes that w(a)= w0 + wa (1-a/a0); but
      Recfast does not assume anything */

  return _SUCCESS_;
}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input/Output: pointer to initialized background structure
 * @return the error status
 */

int background_init(
                    struct precision * ppr,
                    struct background * pba
                    ) {

  /** Summary: */

  /** - define local variables */
  int n_ncdm;
  double rho_ncdm_rel,rho_nu_rel;
  double Neff;
  double w_fld, dw_over_da, integral_fld;
  int filenum=0;

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");

    /* below we want to inform the user about ncdm species*/
    if (pba->N_ncdm > 0) {

      Neff = pba->Omega0_ur/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;

      /* loop over ncdm species */
      for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

        /* inform if p-s-d read in files */
        if (pba->got_files[n_ncdm] == _TRUE_) {
          printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
          filenum++;
        }

        /* call this function to get rho_ncdm */
        background_ncdm_momenta(
                                pba,
                                n_ncdm,
                                pba->q_ncdm_bg[n_ncdm],
                                pba->w_ncdm_bg[n_ncdm],
                                pba->q_size_ncdm_bg[n_ncdm],
                                0.,
                                pba->factor_ncdm[n_ncdm],
                                0.,
                                NULL,
                                &rho_ncdm_rel,
                                NULL,
                                NULL,
                                NULL);

        /* inform user of the contribution of each species to
           radiation density (in relativistic limit): should be
           between 1.01 and 1.02 for each active neutrino species;
           evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
           density of one neutrino in the instantaneous decoupling
           limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
           from the definition of N_eff) */
        rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
          pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4);

        printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
               n_ncdm+1,
               pba->q_size_ncdm_bg[n_ncdm],
               pba->q_size_ncdm[n_ncdm],
               rho_ncdm_rel/rho_nu_rel);

        Neff += rho_ncdm_rel/rho_nu_rel;

      }

      printf(" -> total N_eff = %g (sumed over ultra-relativistic and ncdm species)\n",Neff);

    }
  }

  /** - if shooting failed during input, catch the error here */
  class_test(pba->shooting_failed == _TRUE_,
             pba->error_message,
             "Shooting failed, try optimising input_get_guess(). Error message:\n\n%s",
             pba->shooting_error);

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(pba),
             pba->error_message,
             pba->error_message);

  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  /* Many users asked for this test to be supressed. It is commented out. */
  /*class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
             pba->error_message,
             "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);*/

  class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 -1.)>ppr->smallest_allowed_variation,
             pba->error_message,
             "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",pba->H0,pba->H0/1.e5* _c_,pba->h);

  /* T_cmb in K */
  /* Many users asked for this test to be supressed. It is commented out. */
  /*class_test((pba->T_cmb < _TCMB_SMALL_)||(pba->T_cmb > _TCMB_BIG_),
             pba->error_message,
             "T_cmb=%g out of bounds (%g<T_cmb<%g)",pba->T_cmb,_TCMB_SMALL_,_TCMB_BIG_);*/

  /* Omega_k */
  /* Many users asked for this test to be supressed. It is commented out. */
  /*class_test((pba->Omega0_k < _OMEGAK_SMALL_)||(pba->Omega0_k > _OMEGAK_BIG_),
             pba->error_message,
             "Omegak = %g out of bounds (%g<Omegak<%g) \n",pba->Omega0_k,_OMEGAK_SMALL_,_OMEGAK_BIG_);*/

  /* fluid equation of state */
  if (pba->has_fld == _TRUE_) {

    class_call(background_w_fld(pba,0.,&w_fld,&dw_over_da,&integral_fld), pba->error_message, pba->error_message);

    class_test(w_fld >= 1./3.,
               pba->error_message,
               "Your choice for w(a--->0)=%g is suspicious, since it is bigger than -1/3 there cannot be radiation domination at early times\n",
               w_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
             n_ncdm+1,
             pba->m_ncdm_in_eV[n_ncdm],
             pba->m_ncdm_in_eV[n_ncdm]*pba->deg_ncdm[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);
    }
  }

  /* check other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
             pba->error_message,
             "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
             pba->error_message,
             "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  /** - this function integrates the background over time, allocates
      and fills the background table */
  class_call(background_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  /** - this function finds and stores a few derived parameters at radiation-matter equality */
  class_call(background_find_equality(ppr,pba),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
                    struct background *pba
                    ) {
  int err;

  free(pba->tau_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->background_table);
  free(pba->d2background_dtau2_table);

  err = background_free_input(pba);

  return err;
}

/**
 * Free only the memory space NOT allocated through input_read_parameters()
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free_noinput(
                    struct background *pba
                    ) {
  free(pba->tau_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->background_table);
  free(pba->d2background_dtau2_table);

  /*SI*/
  if ((pba->N_ncdm>0)&&(pba->sum_NR_SI_decoupled>0)) {
      for(int k=0; k<pba->N_ncdm; k++){
      free(pba->ncdm_NR_SI_DF_factor[k]);
      free(pba->ncdm_NR_SI_DF_factor_bg[k]);
      free(pba->dlnf0_dlnq_ncdm_NR_SI[k]);
    }
    free(pba->ncdm_NR_SI_DF_factor);
    free(pba->ncdm_NR_SI_DF_factor_bg);
    free(pba->dlnf0_dlnq_ncdm_NR_SI);
  }
  /*SI*/

  return _SUCCESS_;
}
/**
 * Free pointers inside background structure which were
 * allocated in input_read_parameters()
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_free_input(
                          struct background *pba
                          ) {

  int k;
  if (pba->Omega0_ncdm_tot != 0.){
    for(k=0; k<pba->N_ncdm; k++){
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
    }
    free(pba->ncdm_quadrature_strategy);
    free(pba->ncdm_input_q_size);
    free(pba->ncdm_qmax);
    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->factor_ncdm);
    if(pba->got_files!=NULL)
      free(pba->got_files);
    if(pba->ncdm_psd_files!=NULL)
      free(pba->ncdm_psd_files);
    if(pba->ncdm_psd_parameters!=NULL)
      free(pba->ncdm_psd_parameters);

    /*SI*/
    if(pba->ncdm_si_type!=NULL) free(pba->ncdm_si_type);
    if(pba->ncdm_si_tau_files!=NULL) free(pba->ncdm_si_tau_files);
    
    if(pba->si_ncdm_table_T!=NULL) free(pba->si_ncdm_table_T);
    if(pba->si_ncdm_table_taurel!=NULL) free(pba->si_ncdm_table_taurel);
    if(pba->si_ncdm_table_d2taurel!=NULL) free(pba->si_ncdm_table_d2taurel);
    if(pba->si_ncdm_table_size!=NULL) free(pba->si_ncdm_table_size);

    if(pba->is_early_coupled!=NULL) free(pba->is_early_coupled);
    if(pba->is_NR_SI_decoupled!=NULL) free(pba->is_NR_SI_decoupled);
    if(pba->mass_set_by_DF!=NULL) free(pba->mass_set_by_DF);
    if(pba->Omega_set_by_DF!=NULL) free(pba->Omega_set_by_DF);
    /*SI*/
  }

  if (pba->Omega0_scf != 0.){
    if (pba->scf_parameters != NULL)
      free(pba->scf_parameters);
  }
  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_indices(
                       struct background *pba
                       ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all flags: which species are present? */

  pba->has_cdm = _FALSE_;
  pba->has_ncdm = _FALSE_;
  pba->has_dcdm = _FALSE_;
  pba->has_dr = _FALSE_;
  pba->has_scf = _FALSE_;
  pba->has_lambda = _FALSE_;
  pba->has_fld = _FALSE_;
  pba->has_ur = _FALSE_;
  pba->has_curvature = _FALSE_;
  /*SI*/
  pba->has_si_ncdm = _FALSE_;

  if (pba->Omega0_cdm != 0.)
    pba->has_cdm = _TRUE_;

  if (pba->Omega0_ncdm_tot != 0.)
    pba->has_ncdm = _TRUE_;

  /*SI*/
  // Unused: Sum of the SIWDM type (ie. modes, mode 0 is no SIWDM, mode 1 is interpolation from table)
  int si_type_sum=0;
  for(int i=0; i<pba->N_ncdm; i++) {
    si_type_sum += pba->ncdm_si_type[i];
  }
  if (si_type_sum > 0) pba->has_si_ncdm = _TRUE_;

  if (pba->Omega0_dcdmdr != 0.){
    pba->has_dcdm = _TRUE_;
    if (pba->Gamma_dcdm != 0.)
      pba->has_dr = _TRUE_;
  }

  if (pba->Omega0_scf != 0.)
    pba->has_scf = _TRUE_;

  if (pba->Omega0_lambda != 0.)
    pba->has_lambda = _TRUE_;

  if (pba->Omega0_fld != 0.)
    pba->has_fld = _TRUE_;

  if (pba->Omega0_ur != 0.)
    pba->has_ur = _TRUE_;

  if (pba->sgnK != 0)
    pba->has_curvature = _TRUE_;

  /** - initialize all indices */

  index_bg=0;

  /* index for scale factor */
  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);

  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - index for rho_g (photon density) */
  class_define_index(pba->index_bg_rho_g,_TRUE_,index_bg,1);

  /* - index for rho_b (baryon density) */
  class_define_index(pba->index_bg_rho_b,_TRUE_,index_bg,1);

  /* - index for rho_cdm */
  class_define_index(pba->index_bg_rho_cdm,pba->has_cdm,index_bg,1);

  /* - indices for ncdm. We only define the indices for ncdm1
     (density, pressure, pseudo-pressure), the other ncdm indices
     are contiguous */
  class_define_index(pba->index_bg_rho_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_pseudo_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);

  /* - index for dcdm */
  class_define_index(pba->index_bg_rho_dcdm,pba->has_dcdm,index_bg,1);

  /* - index for dr */
  class_define_index(pba->index_bg_rho_dr,pba->has_dr,index_bg,1);

  /* - indices for scalar field */
  class_define_index(pba->index_bg_phi_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_phi_prime_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_V_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_dV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_ddV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_rho_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_scf,pba->has_scf,index_bg,1);

  /* - index for Lambda */
  class_define_index(pba->index_bg_rho_lambda,pba->has_lambda,index_bg,1);

  /* - index for fluid */
  class_define_index(pba->index_bg_rho_fld,pba->has_fld,index_bg,1);
  class_define_index(pba->index_bg_w_fld,pba->has_fld,index_bg,1);

  /* - index for ultra-relativistic neutrinos/species */
  class_define_index(pba->index_bg_rho_ur,pba->has_ur,index_bg,1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_r,_TRUE_,index_bg,1);

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */

  /**-------------------------------SI BLOCK----------------------------*/

  //Index of Reltime for WDM species 1 in the background vector
  class_define_index(pba->index_bg_taurel_si_ncdm1,pba->has_si_ncdm,index_bg,pba->N_si_ncdm);

  /**-------------------------------SI BLOCK----------------------------*/

  /*    */

  /* - end of indices in the normal vector of background values */
  pba->bg_size_normal = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(pba->index_bg_rho_crit,_TRUE_,index_bg,1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_m,_TRUE_,index_bg,1);

  /* -> conformal distance */
  class_define_index(pba->index_bg_conf_distance,_TRUE_,index_bg,1);

  /* -> angular diameter distance */
  class_define_index(pba->index_bg_ang_distance,_TRUE_,index_bg,1);

  /* -> luminosity distance */
  class_define_index(pba->index_bg_lum_distance,_TRUE_,index_bg,1);

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bg_time,_TRUE_,index_bg,1);

  /* -> conformal sound horizon */
  class_define_index(pba->index_bg_rs,_TRUE_,index_bg,1);

  /* -> density growth factor in dust universe */
  class_define_index(pba->index_bg_D,_TRUE_,index_bg,1);

  /* -> velocity growth factor in dust universe */
  class_define_index(pba->index_bg_f,_TRUE_,index_bg,1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  pba->bg_size = index_bg;

  /* - now, indices in vector of variables to integrate.
     First {B} variables, then {C} variables. */

  index_bi=0;

  /* -> scale factor */
  class_define_index(pba->index_bi_a,_TRUE_,index_bi,1);

  /* -> energy density in DCDM */
  class_define_index(pba->index_bi_rho_dcdm,pba->has_dcdm,index_bi,1);

  /* -> energy density in DR */
  class_define_index(pba->index_bi_rho_dr,pba->has_dr,index_bi,1);

  /* -> energy density in fluid */
  class_define_index(pba->index_bi_rho_fld,pba->has_fld,index_bi,1);

  /* -> scalar field and its derivative wrt conformal time (Zuma) */
  class_define_index(pba->index_bi_phi_scf,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_phi_prime_scf,pba->has_scf,index_bi,1);

  /* End of {B} variables, now continue with {C} variables */
  pba->bi_B_size = index_bi;

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bi_time,_TRUE_,index_bi,1);

  /* -> sound horizon */
  class_define_index(pba->index_bi_rs,_TRUE_,index_bi,1);

  /* -> Second order equation for growth factor */
  class_define_index(pba->index_bi_D,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_D_prime,_TRUE_,index_bi,1);

  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);

  /* -> end of indices in the vector of variables to integrate */
  pba->bi_size = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(pba->index_bi_tau != index_bi-1,
             pba->error_message,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  /* flags for calling the interpolation routine */

  pba->short_info=0;
  pba->normal_info=1;
  pba->long_info=2;

  pba->inter_normal=0;
  pba->inter_closeby=1;

  return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a particular f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx;
  double ksi;
  double qlast,dqlast,f0last,df0last;
  double *param;
  /* Variables corresponding to entries in param: */
  //double square_s12,square_s23,square_s13;
  //double mixing_matrix[3][3];
  //int i;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */
  param = pba->ncdm_psd_parameters; /* extract the optional parameter list from it */
  n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
  ksi = pba->ksi_ncdm[n_ncdm];      /* extract chemical potential */

  /** - shall we interpolate in file, or shall we use analytical formula below? */

  /** - a) deal first with the case of interpolating in files */
  if (pba->got_files[n_ncdm]==_TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];

      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pba->error_message),
                 pba->error_message,     pba->error_message);
    }
  }

  /** - b) deal now with case of reading analytical function */
  else{
    /**
       Next enter your analytic expression(s) for the p.s.d.'s. If
       you need different p.s.d.'s for different species, put each
       p.s.d inside a condition, like for instance: if (n_ncdm==2)
       {*f0=...}.  Remember that n_ncdm = 0 refers to the first
       species.
    */

    /**************************************************/
    /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
    /**************************************************/

    *f0 = 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));

    /**************************************************/

    /** This form is only appropriate for approximate studies, since in
        reality the chemical potentials are associated with flavor
        eigenstates, not mass eigenstates. It is easy to take this into
        account by introducing the mixing angles. In the later part
        (not read by the code) we illustrate how to do this. */

    if (_FALSE_) {

      /* We must use the list of extra parameters read in input, stored in the
         ncdm_psd_parameter list, extracted above from the structure
         and now called param[..] */

      /* check that this list has been read */
      class_test(param == NULL,
                 pba->error_message,
                 "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

      /* extract values from the list (in this example, mixing angles) */
      double square_s12=param[0];
      double square_s23=param[1];
      double square_s13=param[2];

      /* infer mixing matrix */
      double mixing_matrix[3][3];
      int i;

      mixing_matrix[0][0]=pow(fabs(sqrt((1-square_s12)*(1-square_s13))),2);
      mixing_matrix[0][1]=pow(fabs(sqrt(square_s12*(1-square_s13))),2);
      mixing_matrix[0][2]=fabs(square_s13);
      mixing_matrix[1][0]=pow(fabs(sqrt((1-square_s12)*square_s13*square_s23)+sqrt(square_s12*(1-square_s23))),2);
      mixing_matrix[1][1]=pow(fabs(sqrt(square_s12*square_s23*square_s13)-sqrt((1-square_s12)*(1-square_s23))),2);
      mixing_matrix[1][2]=pow(fabs(sqrt(square_s23*(1-square_s13))),2);
      mixing_matrix[2][0]=pow(fabs(sqrt(square_s12*square_s23)-sqrt((1-square_s12)*square_s13*(1-square_s23))),2);
      mixing_matrix[2][1]=pow(sqrt((1-square_s12)*square_s23)+sqrt(square_s12*square_s13*(1-square_s23)),2);
      mixing_matrix[2][2]=pow(fabs(sqrt((1-square_s13)*(1-square_s23))),2);

      /* loop over flavor eigenstates and compute psd of mass eigenstates */
      *f0=0.0;
      for(i=0;i<3;i++){

    	*f0 += mixing_matrix[i][n_ncdm]*1.0/pow(2*_PI_,3)*(1./(exp(q-pba->ksi_ncdm[i])+1.) +1./(exp(q+pba->ksi_ncdm[i])+1.));

      }
    } /* end of region not used, but shown as an example */
  }

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */

int background_ncdm_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * NEW TO SIWDM MODULE
 * 
 * Similar to the function above, this is the function involved in finding
 * optimal quadratures for the case of SOWDM. It is important to modify it
 * to also include for the possible case of a " non relativistic - type"
 * distribution function.
 * 
 * From the documentation above:
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */
int background_ncdm_NR_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  //double c = 0.;
  //double d = 0.;
  //double e = 0.;

  double f = 8./3.*pow(_PI_,0.5) * 3.0;
  double g = (4./pow(_PI_,0.5)) * 3.0;
 
  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q+ f*q*q*q*q*(exp(q-q*q)+exp(-q-q)) + g*q*q*(exp(q-q*q)+exp(-q*q)) );

  return _SUCCESS_;
}

/**
 * SIGNIFICANTLY MODIFIED FOR SIWDM MODULE
 * 
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */
int background_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  /* Allocate pointer arrays: */
  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);

  /* Allocate pointers: */
  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);

  for(k=0, filenum=0; k<pba->N_ncdm; k++){
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
                 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      for (row=0,status=2; status==2; row++){
        status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++){
        status = fscanf(psdfile,"%lf %lf",
                        &pbadist.q[row],&pbadist.f0[row]);
        //		printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);

      //for (int i=0; i<pbadist.tablesize; i++) printf("q= %e, f0=%e, d2f0=%e\n", pbadist.q[i], pbadist.f0[i], pbadist.d2f0[i] );

      filenum++;
    }

    /* Handle perturbation qsampling: */
    if (pba->ncdm_quadrature_strategy[k]==qm_auto){
      /** Automatic q-sampling for this species */
      class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);

      /* If the species self-decouples while non relativistic, 
       * use the new quadrature test function.
       */

      if (pba->is_NR_SI_decoupled[k]){
        class_call(get_qsampling(pba->q_ncdm[k],
  			       pba->w_ncdm[k],
  			       &(pba->q_size_ncdm[k]),
  			       _QUADRATURE_MAX_,
  			       ppr->tol_ncdm,
  			       pbadist.q,
  			       pbadist.tablesize,
  			       background_ncdm_NR_test_function,
  			       background_ncdm_distribution,
  			       &pbadist,
  			       pba->error_message),
  		 pba->error_message,
  		 pba->error_message);
      }
      else{
        class_call(get_qsampling(pba->q_ncdm[k],
               pba->w_ncdm[k],
               &(pba->q_size_ncdm[k]),
               _QUADRATURE_MAX_,
               ppr->tol_ncdm,
               pbadist.q,
               pbadist.tablesize,
               background_ncdm_test_function,
               background_ncdm_distribution,
               &pbadist,
               pba->error_message),
       pba->error_message,
       pba->error_message); 
      }

      pba->q_ncdm[k]=realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
      pba->w_ncdm[k]=realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));

      if (pba->background_verbose > 0){
        if(pba->is_NR_SI_decoupled[k])
        	printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration (taking into account NR SI decoupling)\n",
        	       k+1,
        	       pba->q_size_ncdm[k]);
        else
          printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
                 k+1,
                 pba->q_size_ncdm[k]);
      }

      /* Handle background q_sampling: */
      class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

      /* If the species self-decouples while non relativistic, 
       * use the new quadrature test function.
       */

      if (pba->is_NR_SI_decoupled[k]){
        class_call(get_qsampling(pba->q_ncdm_bg[k],
            pba->w_ncdm_bg[k],
            &(pba->q_size_ncdm_bg[k]),
            _QUADRATURE_MAX_BG_,
            ppr->tol_ncdm_bg,
            pbadist.q,
            pbadist.tablesize,
            background_ncdm_NR_test_function,
            background_ncdm_distribution,
            &pbadist,
            pba->error_message),
        pba->error_message,
        pba->error_message);
      }
      else{
        class_call(get_qsampling(pba->q_ncdm_bg[k],
			        pba->w_ncdm_bg[k],
			        &(pba->q_size_ncdm_bg[k]),
			        _QUADRATURE_MAX_BG_,
			        ppr->tol_ncdm_bg,
			        pbadist.q,
			        pbadist.tablesize,
			        background_ncdm_test_function,
			        background_ncdm_distribution,
			        &pbadist,
			        pba->error_message),
		    pba->error_message,
		    pba->error_message);
      }

      pba->q_ncdm_bg[k]=realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));
      pba->w_ncdm_bg[k]=realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));

      /** - in verbose mode, inform user of number of sampled momenta
	  for background quantities */
      if (pba->background_verbose > 0){
        if (pba->is_NR_SI_decoupled[k])
	        printf("ncdm species i=%d sampled with %d points for purpose of background integration (taking into account NR SI decoupling)\n",
	        k+1,
	        pba->q_size_ncdm_bg[k]);
        else 
          printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
          k+1,
          pba->q_size_ncdm_bg[k]);
      }
    }
    else{
      /** Manual q-sampling for this species. Same sampling used for both perturbation and background sampling, since this will usually be a high precision setting anyway */
      pba->q_size_ncdm_bg[k] = pba->ncdm_input_q_size[k];
      pba->q_size_ncdm[k] = pba->ncdm_input_q_size[k];
      class_alloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double),pba->error_message);
      class_alloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double),pba->error_message);
      class_call(get_qsampling_manual(pba->q_ncdm[k],
				      pba->w_ncdm[k],
				      pba->q_size_ncdm[k],
				      pba->ncdm_qmax[k],
				      pba->ncdm_quadrature_strategy[k],
				      pbadist.q,
				      pbadist.tablesize,
				      background_ncdm_distribution,
				      &pbadist,
				      pba->error_message),
		 pba->error_message,
		 pba->error_message);
      for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
	pba->q_ncdm_bg[k][index_q] = pba->q_ncdm[k][index_q];
	pba->w_ncdm_bg[k][index_q] = pba->w_ncdm[k][index_q];
      }
    /** - in verbose mode, inform user of number of sampled momenta
        for background quantities */
      if (pba->background_verbose > 0)
	printf("ncdm species i=%d sampled with %d points for purpose of background andperturbation integration using the manual method\n",
	       k+1,
	       pba->q_size_ncdm[k]);
    }

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);


    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0),
                 pba->error_message,pba->error_message);

      //Loop to find appropriate dq:
      for(tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){

        if (index_q == 0){
          dq = MIN((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }
        else if (index_q == pba->q_size_ncdm[k]-1){
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

        class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2),
                   pba->error_message,pba->error_message);
        class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2),
                   pba->error_message,pba->error_message);

        if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1),
                 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1),
                 pba->error_message,pba->error_message);
      //5 point estimate of the derivative:
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
      //Avoid underflow in extreme tail:
      if (fabs(f0)==0.)
        pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
        pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
    }

    pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(pba->T_cmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;

    /* If allocated, deallocate interpolation table:  */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }

  /* SIWDM non relativistic DF factors storage:
   * the values of relevant quantities (i.e. f(q) and dlnf/dlnq)
   * are pre-calculated and stored in the deep relativistic and 
   * non-relativistic regimes.
   */

  //printf("Just before NR_SI_block on background_ncdm_init\n");
  if (pba->sum_NR_SI_decoupled>0){
    //printf("Entered NR_SI block on background_ncdm_init\n");
    background_ncdm_NR_SI_DF_store_factors(pba);
  }


  return _SUCCESS_;
}

/**
 * For a given ncdm species: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 * 
 * MODIFIED FOR SIWDM MODULE: An additional normalization
 * factor for the momentum integrals is added in case of a non relativistic
 * DF. This additional factor is either retrieved from memory or calculated,
 * depending on the situation.
 *
 * @param qvec     Input: sampled momenta
 * @param wvec     Input: quadrature weights
 * @param qsize    Input: number of momenta/weights
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redshift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Output: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int background_ncdm_momenta(
                            /* Only calculate for non-NULL pointers: */
                            struct background *pba,
                            int n_ncdm,
                            double * qvec,
                            double * wvec,
                            int qsize,
                            double M,
                            double factor,
                            double z,
                            double * n,
                            double * rho, // density
                            double * p,   // pressure
                            double * drho_dM,  // d rho / d M used in next function
                            double * pseudo_p  // pseudo-p used in ncdm fluid approx
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;
  double factor_prime;
  double T_over_m = ( pba->T_cmb * pba->T_ncdm[n_ncdm] / pba->a_today * ( 1. + z )) * (_k_B_) /
                    ( pba->m_ncdm_in_eV[n_ncdm] *_eV_);
  /** Summary: */

  /** - rescale normalization at given redshift */
  factor2 = factor*pow(1+z,4);

  /** - initialize quantities */
  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  /** - loop over momenta */
  for (index_q=0; index_q<qsize; index_q++) {

    /* squared momentum */
    q2 = qvec[index_q]*qvec[index_q];

    /* energy */
    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    /*Renormalization of DF for NR SI decoupling */
    // Here the additional integration factor is calculated

    //if (pba->is_NR_SI_decoupled[n_ncdm]) background_get_NR_SI_factor_prime(pba,n_ncdm,qvec[index_q],z,&factor_prime);
    //else factor_prime = 1.0;

    //printf("Just before NR_SI block on background_momenta\n");
    if (pba->is_NR_SI_decoupled[n_ncdm]==_TRUE_){

      //printf("Entered NR_SI block on background_momenta\n");
      /*Two different cases for background and perturbation integration...
       *the function background_ncdm_NR_SI_switching switches between ultra realtivistic,
       *non relativistic (precalcualted) and intermediate (calcualted on the fly) regimes
       */
      if (qsize==pba->q_size_ncdm[n_ncdm]){
        //fprintf(stderr, "qsize = %d, T_over_m = %e, pba->ncdm_NR_SI_DF_factor[n_ncdm][index_q] = %e \n", qsize, T_over_m, pba->ncdm_NR_SI_DF_factor[n_ncdm][index_q]);
        background_ncdm_NR_SI_switching(
          T_over_m,
          class_call(background_get_NR_SI_factor_prime(pba,n_ncdm,index_q,z,&factor_prime), pba->error_message, pba->error_message),
          pba->ncdm_NR_SI_DF_factor[n_ncdm][index_q],
          1.0,
          &factor_prime);
      }
      else if (qsize==pba->q_size_ncdm_bg[n_ncdm]){
        background_ncdm_NR_SI_switching(
          T_over_m,
          class_call(background_get_NR_SI_factor_prime_bg(pba,n_ncdm,index_q,z,&factor_prime), pba->error_message, pba->error_message),
          pba->ncdm_NR_SI_DF_factor_bg[n_ncdm][index_q],
          1.0,
          &factor_prime); 
      }
      else{
        class_test(1>0, pba->error_message, "ERROR in q_size in background_ncdm_momenta");
      }

    }
    else{
      factor_prime=1.0;
    }
    //printf("Factor_prime = %e\n", factor_prime );

    /* integrand of the various quantities */
    /* In the case of SIWDM, the additional factor is included, otherwise it's just 1.0*/
    if (n!=NULL) *n += factor_prime*q2*wvec[index_q];
    if (rho!=NULL) *rho += factor_prime*q2*epsilon*wvec[index_q];
    if (p!=NULL) *p += factor_prime*q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM!=NULL) *drho_dM += factor_prime*q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += factor_prime*pow(q2/epsilon,3)/3.0*wvec[index_q];
  }

  /** - adjust normalization */
  if (n!=NULL) *n *= factor2/(1.+z);
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dM!=NULL) *drho_dM *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  return _SUCCESS_;
}

/**
 * When the user passed the density fraction Omega_ncdm or
 * omega_ncdm in input but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

int background_ncdm_M_from_Omega(
                                 struct precision *ppr,
                                 struct background *pba,
                                 int n_ncdm
                                 ) {
  double rho0,rho,n,M,deltaM,drhodM;
  int iter,maxiter=50;

  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[n_ncdm]; /*Remember that rho is defined such that H^2=sum(rho_i) */
  M = 0.0;

  background_ncdm_momenta(
                          pba,
                          n_ncdm,
                          pba->q_ncdm_bg[n_ncdm],
                          pba->w_ncdm_bg[n_ncdm],
                          pba->q_size_ncdm_bg[n_ncdm],
                          M,
                          pba->factor_ncdm[n_ncdm],
                          0.,
                          &n,
                          &rho,
                          NULL,
                          NULL,
                          NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0<rho,pba->error_message,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
             n_ncdm,pba->Omega0_ncdm[n_ncdm],pba->Omega0_ncdm[n_ncdm]*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zeroth order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (iter=1; iter<=maxiter; iter++){

    /* Newton iteration. First get relevant quantities at M: */
    background_ncdm_momenta(
                            pba,
                            n_ncdm,
                            pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            NULL,
                            &rho,
                            NULL,
                            &drhodM,
                            NULL);

    deltaM = (rho0-rho)/drhodM; /* By definition of the derivative */
    if ((M+deltaM)<0.0) deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M)<ppr->tol_M_ncdm){
      /* Accuracy reached.. */
      pba->M_ncdm[n_ncdm] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
             "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

/** ---------------SELF INTERATIONS BLOCK-----------*/
/** These are the functions exclusive to the SIWDM module. 
 * Other existing functions are modified to account for these */

/** 
 * ADDED FOR SIWDM MODULE
 * 
 * Main initializer function for the SIWDM species.
 * This function is only called when (relevant) SIWDM species are detected in the input,
 * and after all initial flags are set. 
 * 
 * For now, only interactions in mode 1 (lookup table) are implemented. For this case, 
 * it is neccesary to allocate the table and retrieve it from file: this function's main purpose
 * is to do that.
 * 
 * The organization of this function is as follows:
 * 
 * - Allocate a few neccesary flags for the SIWDM structures
 * - Regarding species in mode 1:
 * --Allocate relaxation time lookup table
 * --Save it to pba
 * --Test interpolation of the lookup table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */
int background_si_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  if (pba->background_si_verbose>2) {
    printf("%s\n", "Called background_si_ncdm_init");

    printf("ncdm_NR_SI_decoupling_method set to %d \n", pba->ncdm_NR_SI_decoupling_method);
  }

  class_test( ( pba->N_ncdm == 0. ) , pba->error_message, "Self Interacting NCDM called but no NCDM species found.");

  // Allocation and storage of pba->is_early_coupled: this stores whether or not the SIWDM interactions
  // are strong enough to remain coupled at early times, if it isn't the evolution follows WDM (no recoupling!)

  class_alloc(pba->is_early_coupled,sizeof(int)*pba->N_ncdm,pba->error_message);

  /** check if correctly allocated! */

  /*pba->is_early_coupled[pba->N_ncdm-1] = 1;
  if(pba->background_si_verbose>2){
    printf("Checking if correctly allocated: pba->is_early_coupled[pba->N_ncdm] = %d\n", pba->is_early_coupled[pba->N_ncdm-1] );
  }*/

  /* Note from future Rafael: Don't know what this is. Sounds ominous. Afraid to erase it
  FOR NOW, NOT NECCESARY. JUST TESTING. NECESSARY IMPLEMENTATION OF EARLY TCA IN background_initial_conditions(...) */
  /* FOR MODE 1 -> TABLE OF TAU_REL (T_CMB) , INITIALIZE TAU_REL TABLE */

  int k,row,status,n_si_ncdm;
  double tmp1,tmp2;
  double T_test;
  double * taurel_test;
  int tablesize;
  int * last_index;
  FILE *taufile;

  class_alloc(pba->si_ncdm_table_T,sizeof(double*)*pba->N_si_ncdm,pba->error_message);
  class_alloc(pba->si_ncdm_table_taurel,sizeof(double*)*pba->N_si_ncdm,pba->error_message);
  class_alloc(pba->si_ncdm_table_d2taurel,sizeof(double*)*pba->N_si_ncdm,pba->error_message);
  class_alloc(pba->si_ncdm_table_size,sizeof(double)*pba->N_si_ncdm,pba->error_message);

  for(k=0, n_si_ncdm=0; k<pba->N_ncdm; k++){

    /* CASE: MODE 1. DO WE HAVE A TABLE TO INTERPOLATE? */
    if (pba->ncdm_si_type[k]==1){

      if(pba->background_si_verbose>1){
        printf("NCDM species %d found of interaction type 1.\n", k+1);
      }

      //try to open file
      taufile = fopen(pba->ncdm_si_tau_files+n_si_ncdm*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(taufile == NULL,pba->error_message,
                 "Could not open file %s for ncdm species %d!",pba->ncdm_si_tau_files+n_si_ncdm*_ARGUMENT_LENGTH_MAX_, k);

      if(pba->background_si_verbose>1){
        printf("Successfuly opened file %s\n", pba->ncdm_si_tau_files+n_si_ncdm*_ARGUMENT_LENGTH_MAX_);
      }

      // Find size of table:
      for (row=0,status=2; status==2; row++){
        status = fscanf(taufile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(taufile);
      tablesize = row-1;

      if(pba->background_si_verbose>1){
        printf("Table size found to be %d. Allocating memory.\n", tablesize );
      }

      /*Allocate room for interpolation table: */
      class_alloc(pba->si_ncdm_table_T[n_si_ncdm],sizeof(double)*tablesize,pba->error_message);
      class_alloc(pba->si_ncdm_table_taurel[n_si_ncdm],sizeof(double)*tablesize,pba->error_message);
      class_alloc(pba->si_ncdm_table_d2taurel[n_si_ncdm],sizeof(double)*tablesize,pba->error_message);

      if(pba->background_si_verbose>1){
        printf("Memory allocated. Writing...\n");
      }

      //Store Table
      for (row=0; row<tablesize; row++){
        status = fscanf(taufile,"%lf %lf",
                        &(pba->si_ncdm_table_T[n_si_ncdm][row]),&(pba->si_ncdm_table_taurel[n_si_ncdm][row]));
        if(pba->background_si_verbose>3){
          printf("(T,Tau_rel) = (%g,%g)\n",pba->si_ncdm_table_T[n_si_ncdm][row],pba->si_ncdm_table_taurel[n_si_ncdm][row]);
        }
      }
      fclose(taufile);

      if(pba->background_si_verbose>1){
        printf("Stored Table. T[last] = %g, Tau_rel[last] = %g \n", pba->si_ncdm_table_T[n_si_ncdm][tablesize-1], pba->si_ncdm_table_taurel[n_si_ncdm][tablesize-1]);
      }

      //printf("Memory allocated. Writing...\n");

      /* It is much better to store and interpolate from the log() of these quantities.
       * We do that for accuracy of the interpolation 
       */

      /* Setting up the log() of T, tau_rel in the tables. Remember to do exp() when evaluating the spline! */

      double temp;

      // Change everything to log()
      for (row=0; row<tablesize; row++){
        temp = pba->si_ncdm_table_T[n_si_ncdm][row];
        pba->si_ncdm_table_T[n_si_ncdm][row] = log(temp);
        temp = pba->si_ncdm_table_taurel[n_si_ncdm][row];
        pba->si_ncdm_table_taurel[n_si_ncdm][row] = log(temp);
        if(pba->background_si_verbose>3){
          printf("(log T,log Tau_rel) = (%g,%g)\n",pba->si_ncdm_table_T[n_si_ncdm][row],pba->si_ncdm_table_taurel[n_si_ncdm][row]);
        }
      }

      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pba->si_ncdm_table_T[n_si_ncdm],
                                          tablesize,
                                          pba->si_ncdm_table_taurel[n_si_ncdm],
                                          1,
                                          pba->si_ncdm_table_d2taurel[n_si_ncdm],
                                          _SPLINE_NATURAL_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);

      // Testing: Output to console

      if(pba->background_si_verbose>2){
        printf("Spline interpolation called.\n");
      }

      if(pba->background_si_verbose>4){
        printf("Printing interpolation table (again) \n");
      }

      for(int i=0; i<tablesize; i++){
        if(pba->background_si_verbose>4){
          printf("(log T, log taurel, log d2taurel) = (%e,%e,%e)\n", pba->si_ncdm_table_T[n_si_ncdm][i], pba->si_ncdm_table_taurel[n_si_ncdm][i], pba->si_ncdm_table_d2taurel[n_si_ncdm][i] );
        }
      }

      //Store size of table
      pba->si_ncdm_table_size[n_si_ncdm] = tablesize;

      /* In case of any problems, we test the interpolation here.
       * Simply print out to console table values vs interpolation at sampled T
       */

      if(pba->background_si_verbose>1){
        printf("Testing interpolation...\n");
      }

      class_alloc(taurel_test, sizeof(double), pba->error_message);
      class_alloc(last_index, sizeof(int), pba->error_message);
      *taurel_test = 0.0;
      *last_index = 0;

      if(pba->background_si_verbose>4){
        printf("Assigned!\n");
      }

      for (int i=1; i<(tablesize-1); i++){

        T_test = ( pba->si_ncdm_table_T[n_si_ncdm][i-1] + pba->si_ncdm_table_T[n_si_ncdm][i+1] )/2. ;

        class_call(array_interpolate_spline(
          pba->si_ncdm_table_T[n_si_ncdm],
          pba->si_ncdm_table_size[n_si_ncdm],
          pba->si_ncdm_table_taurel[n_si_ncdm],
          pba->si_ncdm_table_d2taurel[n_si_ncdm],
          1,
          T_test,
          last_index,
          taurel_test,
          1,
          pba->error_message),
        pba->error_message,     pba->error_message);

        if(pba->background_si_verbose>3) {
          printf("Test interpolation results: T_test = %g, taurel_test = %g, last_index = %d (interpolated) \n", 
                exp(T_test), exp(*taurel_test), *last_index);
        }
      }

      n_si_ncdm++;
    }

    else{
      class_test( (pba->ncdm_si_type[k]>=2) , pba->error_message, "Self Interaction type 2, 3 not yet implemented!");
    }

  }

  return _SUCCESS_;
}

/** 
 * ADDED FOR SIWDM MODULE
 * 
 * Calculate the relaxation time for ncdm species n_ncdm, at redshift z.
 * If the interactions are set to mode 1: perform spline interpolation on the tau_rel(T_gamma) table stored in pba.
 *
 * @param pba Input: precision structure
 * @param reltime Output: relaxation time
 * @param z Input: Redshift
 * @param n_ncdm Input: si ncdm species number 
 */
int background_si_ncdm_reltime(
                            struct background *pba,
                            double * reltime,
                            double z,
                            int n_ncdm){

  //double T_at_z, T_last, reltime_last, dT_last, dreltime_last;
  double T_at_z, logT_last, logReltime_last, d_logT, d_logReltime, log_slope;
  int last_index, n_si_ncdm;
  double temp_reltime;

  T_at_z = pba->T_cmb * pba->T_ncdm[n_si_ncdm] * (1+z);
  //printf("%e %e\n",z, T_at_z );
  last_index = pba->si_ncdm_table_size[n_si_ncdm] - 1;

  //If the species is SIWDM of mode 1:
  if (pba->ncdm_si_type[n_ncdm]==1){

    //(Important! Different from n_ncdm! Mapping between the two is stored in pba->ncdm_si_index) 
    n_si_ncdm = pba->ncdm_si_index[n_ncdm];

    //Handle T > T_max: tau_rel set to the first T value of the table (ie the highest temperature)
    if(T_at_z > exp( pba->si_ncdm_table_T[n_si_ncdm][0]) ) {
      //Handle T > T_max (first) case :
      //printf("T>T_max. T_at_z = %e, T[0]= %e\n", T_at_z, pba->si_ncdm_table_T[n_si_ncdm][0]);
      *reltime = exp( pba->si_ncdm_table_taurel[n_si_ncdm][0] );
    }
    else if(T_at_z < exp( pba->si_ncdm_table_T[n_si_ncdm][last_index]) ) {
      //printf("T>T_max. last_index= %d, T_at_z = %e, T[last]= %e\n", last_index, T_at_z, pba->si_ncdm_table_T[n_ncdm][last_index]);
      //Handle T < Tmin case : Fit power law to last table values, get value from fit

      logT_last = pba->si_ncdm_table_T[n_si_ncdm][last_index] ;
      logReltime_last =  pba->si_ncdm_table_taurel[n_si_ncdm][last_index] ;
      d_logT= logT_last - pba->si_ncdm_table_T[n_si_ncdm][last_index-1] ;
      d_logReltime = logReltime_last - pba->si_ncdm_table_taurel[n_si_ncdm][last_index-1] ;

      log_slope = d_logReltime / d_logT;

      temp_reltime = exp( ( log(T_at_z) - logT_last ) * log_slope + logReltime_last );

      //If the fit gives a weird value (ie inf, nan) just use the last table value
      if ( isnormal(temp_reltime) ){
        *reltime = temp_reltime;  //Power law approximation
      }
      else{
        *reltime = exp( pba->si_ncdm_table_taurel[n_si_ncdm][last_index] );
      }
    }
    else{
      //Otherwise, do the regular interpolation
      //printf("Interpolating...\n");
      class_call(array_interpolate_spline(
        pba->si_ncdm_table_T[n_si_ncdm],
        pba->si_ncdm_table_size[n_si_ncdm],
        pba->si_ncdm_table_taurel[n_si_ncdm],
        pba->si_ncdm_table_d2taurel[n_si_ncdm],
        1,
        log(T_at_z),
        &last_index,
        &temp_reltime,
        1,
        pba->error_message),
      pba->error_message,     pba->error_message);

      //printf("T[0] = %e, taurel[0] = %e, log(T_at_z)= %e \n", 
      //  exp(pba->si_ncdm_table_T[n_si_ncdm][0]), exp(pba->si_ncdm_table_taurel[n_si_ncdm][0]), log(T_at_z));

      //printf("T_at_z=%e, temp_reltime=%e \n", T_at_z, temp_reltime);

      *reltime = exp(temp_reltime);

    }

  }

  return _SUCCESS_;
}

/** 
 * ADDED FOR SIWDM MODULE
 * 
 * Calculates the temperature and normalization of an equivalent Fermi-Dirac distribution at early times,
 * in case a non-thermal form is given.
 * 
 * This function is called in the background_initial_conditions block. This is done in the case a 
 * background distribution function is given that is not equilibrium (in this case we only consider FD),
 * but we detect that indeed self interactions are relevant at the earliest computed time.
 * A problem arises because, while we have given a non thermal DF, self intearctions should very quickly thermalize the 
 * whole distribution, ie bring it to an equiibrium form. As the self interactions do not change total number density,
 * energy density or any collisional invariants, we construct a new FD distribution that keeps these intact.
 * 
 * The final distribution is simply a Fermi Dirac distribution, characterized by a new temperature and overall 
 * normalization factor, such as f0(q) =  equivalent_factor * (previous normalization) * 1/(exp(q/equivalent_T_ncdm_0)+1).
 * The routine outputs the values for these equivalent quantities.
 *
 * @param ppr Input: precision structure
 * @param pba Input: background structure
 * @param equivalent_T_ncdm_0 Output: Equivalent Fermi-Dirac temperature (in units of T_cmb)
 * @param equivalent_factor Output: Equivalent Fermi-Dirac normalization (in units of the previous normalization factor for f_FD)
 * @param z Input: Redshift
 * @param n_ncdm Input: ncdm species number 
 */ 

int background_si_ncdm_equivalent_params(
                                      struct precision *ppr,
                                      struct background *pba,
                                      double * equivalent_T_ncdm_0,
                                      double * equivalent_factor,
                                      double z,
                                      int n_ncdm){


  double rho;
  double n;

  class_call(background_ncdm_momenta(
   pba,
   n_ncdm,
   pba->q_ncdm_bg[n_ncdm],
   pba->w_ncdm_bg[n_ncdm],
   pba->q_size_ncdm_bg[n_ncdm],
   pba->M_ncdm[n_ncdm],
   pba->factor_ncdm[n_ncdm],
   z,
   &n,
   &rho,
   NULL,
   NULL,
   NULL),
  pba->error_message,
  pba->error_message);

  //printf("Called background_ncdm_momenta...\n");

  n *= 1.0 / (pba->T_cmb * pba->T_ncdm[n_ncdm] * _k_B_) ; // Correcting a bug in the calculation of n_ncdm 

  //printf("rho= %g, n=%g \n", rho, n);

  double k_bT;
  double C;
  double K_3 = 5.6822; // F_3(z=0) the complete fermi dirac integral of order 3 with 0 fugacity
  double K_2 = 1.8030; // F_2(z=0) the complete fermi dirac integral of order 2 with 0 fugacity

  //printf("Defined parameters k_bT, C, K_3, K_2...\n");

  k_bT = rho / n * K_2 / K_3;
  C = rho / pow( k_bT, 4) / K_3;

  //printf("Calculated k_bT=%g, C= %g \n", k_bT, C);

  if (equivalent_T_ncdm_0!=NULL) *equivalent_T_ncdm_0 = k_bT / _k_B_ / pba->T_cmb / (1+z); 
  if (equivalent_factor!=NULL) 
    *equivalent_factor = C / pba->factor_ncdm[n_ncdm] * pow( pba->T_cmb * pba->T_ncdm[n_ncdm] * _k_B_ , 4 ) ;

  if (pba->background_si_verbose>0){
    printf("Finished background_si_ncdm_equivalent_params. New parameters for Fermi_dirac equivalent DF found to be: k_bT=%e C=%e\n",
      k_bT,
      C);
  }

 return _SUCCESS_;
}

/* ADDED FOR SIWDM MODULE 
 * 
 * This function calculates the correction factors in the integration rules for the case of SIWDM,
 * particularly in the case of non relativistic self decoupling.
 *
 * In this case, it is neccesary to correct the integration rules of the form $\int g f_0 \sim \sum g(q_i) w_i
 * in order to account for the change in the background distribution f_0 with time! We achieve
 * this by changing these weights through a correction factor $\sim f_0(q_i,t)/f_0(q_i,t_init)$.
 *
 * In this case, we use as the "true" background distribution a switching function between the 
 * relativistic and nonrelativisitic regimes. This function is specially designed for the perturabtions
 * integration rule.
 *
 * @param pba Input: background structure
 * @param n_ncdm Input: ncdm species number 
 * @param index_q Input: index corresponding to q_i
 * @param z Input: redshift
 * @param factor_prime Output: f_0 correction factor for w_i at redshift z
 */

int background_get_NR_SI_factor_prime(
                                      struct background *pba,
                                      int n_ncdm,
                                      int index_q,
                                      double z,
                                      double *factor_prime
                                      ){

  double T_over_m = ( pba->T_cmb * pba->T_ncdm[n_ncdm] / pba->a_today * ( 1. + z )) * (_k_B_) /
                    ( pba->m_ncdm_in_eV[n_ncdm] *_eV_);
  double q = pba->q_ncdm[n_ncdm][index_q];
  //double switch_f = 

  //If we are in the deep relativistic regime  
  if (T_over_m > _NCDM_LARGE_T_OVER_M_){
    *factor_prime = 1.0;
  }
  else{
    //If we are not, check the options for decoupling method and output
    if ((pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK)||(pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy)){
      //double a = 1 + 3.534 * 0.5 * erfc( log( T_over_m ) );
      //double b = 1 + 0.0748 * 0.5 * erfc( log( T_over_m ) );
      //double c = 1 + 0.5 * erfc( log( T_over_m ) );

      double rel_df = 2.0/pow(2*_PI_,3)*(1./(exp(q)+1.));
      double nr_df;
      if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK){ nr_df = 4.534 * exp(-1.0748*pow(q,2)); }
      if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy){ nr_df = 2.199 * exp(-1.5072*pow(q,2)); }

      //*factor_prime = a*(exp(q-b*pow(q,c))+exp(-b*pow(q,c))) ;
      *factor_prime = 0.5*erfc( log(T_over_m) )*nr_df/rel_df + ( 1.0 - 0.5*erfc( log(T_over_m) ) ) * 1.0 ;
    }
    if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_ignore){
      *factor_prime=1.0;
    }
  }

  //fprintf(stderr, "background_get_NR_SI_factor_prime called for ncdm species %d, momenta %e at t/m = %e. It will output a value of %e, compared to a rel value of %e. \n",
  //n_ncdm, q, T_over_m, *factor_prime, 1.0);
  
  //*factor_prime = 1.0;
  return _SUCCESS_;

}

/* ADDED FOR SIWDM MODULE 
 * 
 * This function calculates the correction factors in the integration rules for the case of SIWDM,
 * particularly in the case of non relativistic self decoupling.
 *
 * In this case, it is neccesary to correct the integration rules of the form $\int g f_0 \sim \sum g(q_i) w_i
 * in order to account for the change in the background distribution f_0 with time! We achieve
 * this by changing these weights through a correction factor $\sim f_0(q_i,t)/f_0(q_i,t_init)$.
 *
 * In this case, we use as the "true" background distribution a switching function between the 
 * relativistic and nonrelativisitic regimes. This function is specially designed for the background
 * integration rule.
 *
 * @param pba Input: background structure
 * @param n_ncdm Input: ncdm species number 
 * @param index_q Input: index corresponding to q_i
 * @param z Input: redshift
 * @param factor_prime Output: f_0 correction factor for w_i at redshift z
 */

int background_get_NR_SI_factor_prime_bg(
                                      struct background *pba,
                                      int n_ncdm,
                                      int index_q,
                                      double z,
                                      double *factor_prime
                                      ){

  double T_over_m = ( pba->T_cmb * pba->T_ncdm[n_ncdm] / pba->a_today * ( 1. + z )) * (_k_B_) /
                    ( pba->m_ncdm_in_eV[n_ncdm] *_eV_);
  double q = pba->q_ncdm_bg[n_ncdm][index_q];

  //If we are in the deep relativistic regime
  if (T_over_m > _NCDM_LARGE_T_OVER_M_){
    *factor_prime = 1.0;
    //*factor_prime = (1.0 + exp(-q))  * (pow(2*_PI_,3) / 2.0) * (7.0 * pow(_PI_,4) / 720. );
  }
  else{
    //If we are not, check the options for decoupling method and output
    if ((pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK)||(pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy)){
      //double a = 1 + 3.534 * 0.5 * erfc( log( T_over_m ) );
      //double b = 1 + 0.0748 * 0.5 * erfc( log( T_over_m ) );
      //double c = 1 + 0.5 * erfc( log( T_over_m ) );

      double nr_df;
      double rel_df = 2.0/pow(2*_PI_,3)*(1./(exp(q)+1.));

      if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK){ nr_df = 4.534 * exp(-1.0748*pow(q,2)); }
      if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy){ nr_df = 2.199 * exp(-1.5072*pow(q,2)); }

      //*factor_prime = a*(exp(q-b*pow(q,c))+exp(-b*pow(q,c))) ;
      *factor_prime = 0.5*erfc( log(T_over_m) )*nr_df/rel_df + ( 1.0 - 0.5*erfc( log(T_over_m) ) ) * 1.0 ;
    }
    if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_ignore){
      *factor_prime=1.0;
    }
  }

  //fprintf(stderr, "background_get_NR_SI_factor_prime_bg called for ncdm species %d, momenta %e at T/m = %e. It will output a value of %e, compared to a rel value of %e. \n",
  //n_ncdm, q, T_over_m, *factor_prime, 1.0);

  return _SUCCESS_;

}

/* ADDED FOR SIWDM MODULE 
 * 
 * This function calculates d ln f_0 / d ln q at q_i in the case of non relativistic self decoupling. 
 * As was the case with the integration rule, this is also a quantity that used to be static but 
 * now changes with time, so this routine calculates these in all regimes.
 *
 * In this case, we use as the "true" background distribution a switching function between the 
 * relativistic and nonrelativisitic regimes. This function is specially designed for the perturbations.
 *
 * @param pba Input: background structure
 * @param n_ncdm Input: ncdm species number 
 * @param index_q Input: index corresponding to q_i
 * @param z Input: redshift
 * @param dlnf0_dlnq Output
 * 
 */

int background_get_NR_SI_dlnf0_dlnq(
                                    struct background *pba,
                                    int n_ncdm,
                                    int index_q,
                                    double z,
                                    double * dlnf0_dlnq
                                    ){

  double T_over_m = ( pba->T_cmb * pba->T_ncdm[n_ncdm] / pba->a_today * ( 1. + z )) * (_k_B_) /
                    ( pba->m_ncdm_in_eV[n_ncdm] *_eV_);
  double q = pba->q_ncdm[n_ncdm][index_q];

  //Only do this if the species is si, actually decouples NR and is not ignored (sanity checks)

  if (pba->ncdm_si_type[n_ncdm] >= 1 && pba->is_NR_SI_decoupled[n_ncdm] == _TRUE_ && pba->ncdm_NR_SI_decoupling_method != (int) NR_SI_ignore){

    // If in the deep relativistic regime, same as before
    if (T_over_m > _NCDM_LARGE_T_OVER_M_){
      //*dlnf0_dlnq = -q;
      *dlnf0_dlnq= - q/(1.0+exp(-q));
    }
    else{
      // Otherwise, calculate
      if ((pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK)||(pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy)){
        //double a = 1 + 3.534 / 2. * erfc( log( T_over_m ) );
        //double b = 1 + 0.0748 * 0.5 * erfc( log( T_over_m ) );
        //double c = 1 + 0.5 * erfc( log( T_over_m ) );

        double nr_df;
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK){ nr_df = 4.534 * exp(-1.0748*pow(q,2)); }
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy){ nr_df = 2.199 * exp(-1.5072*pow(q,2)); }
        //printf("nrdf=%e\n",nr_df);

        double rel_df = 2.0/pow(2*_PI_,3)*(1./(exp(q)+1.));
        //printf("rel_df=%e\n",rel_df);

        double f0 = 0.5*erfc( log(T_over_m) )*nr_df + (1.0-0.5*erfc( log(T_over_m) ) ) *rel_df;
        //printf("f0=%e\n",f0);

        //First calculate df0/dq
        double df0dq = - ( 0.5*erfc( log(T_over_m) )*nr_df*(2.1496*q) + (1.0 - 0.5*erfc( log(T_over_m) ) )*rel_df*(1.0/(1.0+exp(-q))) );
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK){ 
          df0dq = - ( 0.5*erfc( log(T_over_m) )*nr_df*(2.1496*q) + (1.0 - 0.5*erfc( log(T_over_m) ) )*rel_df*(1.0/(1.0+exp(-q))) );
        }
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy){ 
          df0dq = - ( 0.5*erfc( log(T_over_m) )*nr_df*(2*1.5072*q) + (1.0 - 0.5*erfc( log(T_over_m) ) )*rel_df*(1.0/(1.0+exp(-q))) );
        }
        //printf("df0dq=%e\n",df0dq);

        //*dlnf0_dlnq = - b * c * pow(q,c);

        //We may find that some of these values (especially f_0(q)) can be zero. If they are not, business as usual.
        if ((q!=0.0)&&(f0!=0.0)&&(df0dq!=0.0)){
          *dlnf0_dlnq = (q/f0) * df0dq;
        }
        else{
          if (T_over_m<1.0){
            //Assume that we are on the nr case...
            if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_microK){ *dlnf0_dlnq= - 2.0 * 1.0748 * pow(q,2); }
            if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_extrapolateDF_iEntropy){ *dlnf0_dlnq= - 2.0 * 1.5072 * pow(q,2); }
          }
          else if (T_over_m>1.0){
            //Assume that the distribution is FD
            *dlnf0_dlnq = - q / (1 + exp(-q));
          }
        }
      }
    }
  }
  else{
    *dlnf0_dlnq = pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
  }

  //fprintf(stderr, "background_get_NR_SI_dlnf0_dlnq called for ncdm species %d, momenta %e at T/m = %e. It will output a value of %e, compared to a rel value of %e. \n",
  //  n_ncdm, q, T_over_m, *dlnf0_dlnq, pba->dlnf0_dlnq_ncdm[n_ncdm][index_q]);

  return _SUCCESS_;
}

/* ADDED FOR SIWDM MODULE 
 * 
 * This function allocates memory and stores all of these previous quantities (integration correction
 * factors, dlnf0_dlnq) into the background structure in the non realtivistic limitng case, in the case
 * of non relativistic self decoupling. This saves up some computing time when calculating 
 * perturbation evolution.
 *
 * @param pba Input/Output: background structure
 */

int background_ncdm_NR_SI_DF_store_factors(
                                           struct background *pba
                                           ){

  class_test((pba->ncdm_NR_SI_decoupling_method < (int) NR_SI_ignore)||(pba->ncdm_NR_SI_decoupling_method > (int) NR_SI_extrapolateDF_iEntropy),
             pba->error_message,
             "Invalid NR_SI extrapolate method!");

  //Allocate space for these quantities
  class_alloc(pba->ncdm_NR_SI_DF_factor, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->ncdm_NR_SI_DF_factor_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm_NR_SI, sizeof(double*)*pba->N_ncdm,pba->error_message);

  for (int n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){
    class_alloc(pba->ncdm_NR_SI_DF_factor[n_ncdm], pba->q_size_ncdm[n_ncdm]*sizeof(double),pba->error_message);
    class_alloc(pba->ncdm_NR_SI_DF_factor_bg[n_ncdm], pba->q_size_ncdm_bg[n_ncdm]*sizeof(double),pba->error_message);
    class_alloc(pba->dlnf0_dlnq_ncdm_NR_SI[n_ncdm], pba->q_size_ncdm[n_ncdm]*sizeof(double),pba->error_message);

    //Perturbations: we store factor_prime and dlnf0_dlnq
    for (int index_q=0; index_q<pba->q_size_ncdm[n_ncdm]; index_q++){

      // Only store if it is in nr self decoupling case
      if (pba->is_NR_SI_decoupled[n_ncdm]==_TRUE_){       
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_ignore){
          pba->ncdm_NR_SI_DF_factor[n_ncdm][index_q]=1.0;
        }
        else{
          class_call(background_get_NR_SI_factor_prime(pba, n_ncdm, index_q, 0.0, &(pba->ncdm_NR_SI_DF_factor[n_ncdm][index_q])), pba->error_message, pba->error_message);
        }
      }
      else{
        pba->ncdm_NR_SI_DF_factor[n_ncdm][index_q]=1.0;
      }

      if (pba->is_NR_SI_decoupled[n_ncdm]==_TRUE_){
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_ignore){
          pba->dlnf0_dlnq_ncdm_NR_SI[n_ncdm][index_q]=pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
        }
        else{
          class_call(background_get_NR_SI_dlnf0_dlnq(pba, n_ncdm, index_q, 0.0, &(pba->dlnf0_dlnq_ncdm_NR_SI[n_ncdm][index_q])), pba->error_message, pba->error_message);
          if (pba->background_si_verbose>3){
            printf("Inside background_ncdm_NR_SI_DF_store_factors: dln... = %e\n", pba->dlnf0_dlnq_ncdm_NR_SI[n_ncdm][index_q]);
            printf("Just to compare: for the rel case %e\n", pba->dlnf0_dlnq_ncdm[n_ncdm][index_q]);
          }
        }
      }
      else{
        pba->dlnf0_dlnq_ncdm_NR_SI[n_ncdm][index_q]=pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
      }

    }

    //Background: we store only factor_prime
    for (int index_q_bg=0; index_q_bg<pba->q_size_ncdm_bg[n_ncdm]; index_q_bg++){

      // Only store if it is in nr self decoupling case
      if (pba->is_NR_SI_decoupled[n_ncdm]==_TRUE_){ 
        if (pba->ncdm_NR_SI_decoupling_method == (int) NR_SI_ignore){
          pba->ncdm_NR_SI_DF_factor_bg[n_ncdm][index_q_bg]=1.0;
        }
        else{
          class_call(background_get_NR_SI_factor_prime_bg(pba, n_ncdm, index_q_bg, 0.0, &(pba->ncdm_NR_SI_DF_factor_bg[n_ncdm][index_q_bg])), pba->error_message, pba->error_message); 
        }
      }
      else{
        pba->ncdm_NR_SI_DF_factor_bg[n_ncdm][index_q_bg]=1.0;
      }
    }

  }

  return _SUCCESS_;
}

/**------------END OF SELF INTERACTIONS BLOCK-------*/


/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_solve(
                     struct precision *ppr,
                     struct background *pba
                     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw;
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double tau_start;
  /* final conformal time */
  double tau_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  /* vector of all background quantities */
  double * pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index=0;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius=0.;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since tau is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);

  /* here tau_end is in fact the initial time (in the next loop
     tau_start = tau_end) */
  tau_end=pvecback_integration[pba->index_bi_tau];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);

  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps: call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of tau */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    //printf("Called background_functions at background_solve (calculations)\n");
    class_call(background_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pba->error_message);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
                               pData,
                               pba->bi_size,
                               pba->bt_size,
                               pba->index_bi_a,
                               pba->a_today,
                               &last_index,
                               pvecback_integration,
                               pba->bi_size,
                               pba->error_message),
             pba->error_message,
             pba->error_message);

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++)
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  if (pba->has_dcdm == _TRUE_){
    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_){
    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
  }


  /** - allocate background tables */
  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2tau_dz2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_dtau2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */
  for (i=0; i < pba->bt_size; i++) {

    /* -> establish correspondence between the integrated variable and the bg variables */

    pba->tau_table[i] = pData[i*pba->bi_size+pba->index_bi_tau];

    class_test(pData[i*pba->bi_size+pba->index_bi_a] <= 0.,
               pba->error_message,
               "a = %e instead of strictly positiv",pData[i*pba->bi_size+pba->index_bi_a]);

    pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;

    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_tau];

    if (pba->sgnK == 0) comoving_radius = pvecback[pba->index_bg_conf_distance];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(-pba->K);

    pvecback[pba->index_bg_ang_distance] = pba->a_today*comoving_radius/(1.+pba->z_table[i]);
    pvecback[pba->index_bg_lum_distance] = pba->a_today*comoving_radius*(1.+pba->z_table[i]);
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities depending only on {B} variables.
       The value of {B} variables in pData are also copied to pvecback.*/
    //printf("Called background_functions at background_solve (store in table)\n");
    class_call(background_functions(pba,pData+i*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

    /* -> compute growth functions (valid in dust universe) */

    /* Normalise D(z=0)=1 and construct f = D_prime/(aHD) */
    pvecback[pba->index_bg_D] = pData[i*pba->bi_size+pba->index_bi_D]/pData[(pba->bt_size-1)*pba->bi_size+pba->index_bi_D];
    pvecback[pba->index_bg_f] = pData[i*pba->bi_size+pba->index_bi_D_prime]/
      (pData[i*pba->bi_size+pba->index_bi_D]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);

    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_table");
  }

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);

  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      1,
                                      pba->d2tau_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table,
                                      pba->bt_size,
                                      pba->background_table,
                                      pba->bg_size,
                                      pba->d2background_dtau2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  /** - compute remaining "related parameters"
   *     - so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */
  pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
               *pba->background_table[pba->index_bg_rho_crit]
               -pba->background_table[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pba->background_table[pba->index_bg_rho_g]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);
  }

  if (pba->background_verbose > 2) {
    if ((pba->has_dcdm == _TRUE_)&&(pba->has_dr == _TRUE_)){
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n",pba->Omega0_dcdm);
      printf("     -> Omega0_dr = %f\n",pba->Omega0_dr);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
             pba->Omega0_dr+pba->Omega0_dcdm,pba->Omega0_dcdmdr);
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n",pba->Omega_ini_dcdm/pba->Omega0_b);
    }
    if (pba->has_scf == _TRUE_){
      printf("    Scalar field details:\n");
      printf("     -> Omega_scf = %g, wished %g\n",
             pvecback[pba->index_bg_rho_scf]/pvecback[pba->index_bg_rho_crit], pba->Omega0_scf);
      if(pba->has_lambda == _TRUE_)
	printf("     -> Omega_Lambda = %g, wished %g\n",
               pvecback[pba->index_bg_rho_lambda]/pvecback[pba->index_bg_rho_crit], pba->Omega0_lambda);
      printf("     -> parameters: [lambda, alpha, A, B] = \n");
      printf("                    [");
      for (i=0; i<pba->scf_parameters_size-1; i++){
        printf("%.3f, ",pba->scf_parameters[i]);
      }
      printf("%.3f]\n",pba->scf_parameters[pba->scf_parameters_size-1]);
    }
  }

  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 * 
 * MODIFIED FOR SIWDM MODULE
 * At this function (before integration) we perform many of the initial checks 
 * that are needed for internal consistency of the SIWDM module.
 * 
 * Why do we do these checks here? 
 * To correctly account for SIWDM, we need to determine if the species is coupled
 * at the earliest time, if it decouples while non relativistic, etc.
 * It is hard to perform these checks without the (non cold)DM variables already 
 * defined, so the approach taken is to first define the potential SIWDM species
 * as standard WDM. The program will attempt to begin integration and reach the 
 * SIWDM block in this function, where we then check what conditions about thermal 
 * history does it fulfill, re-define the cosmological species according 
 * to this thermal history and keep going.
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities used as workspace
 * @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
 * @return the error status
 */

int background_initial_conditions(
                                  struct precision *ppr,
                                  struct background *pba,
                                  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
                                  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  double rho_ncdm, p_ncdm, rho_ncdm_rel_tot=0.;
  double f,Omega_rad, rho_rad;
  int counter,is_early_enough,n_ncdm;
  double scf_lambda;
  double rho_fld_today;
  double w_fld,dw_over_da_fld,integral_fld;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today;

  /**  If we have ncdm species, perhaps we need to start earlier
      than the standard value for the species to be relativistic.
      This could happen for some WDM models.
  */

  if (pba->has_ncdm == _TRUE_) {

    for (counter=0; counter < _MAX_IT_; counter++) {

      is_early_enough = _TRUE_;
      rho_ncdm_rel_tot = 0.;

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

	class_call(background_ncdm_momenta(
             pba,
             n_ncdm,
             pba->q_ncdm_bg[n_ncdm],
					   pba->w_ncdm_bg[n_ncdm],
					   pba->q_size_ncdm_bg[n_ncdm],
					   pba->M_ncdm[n_ncdm],
					   pba->factor_ncdm[n_ncdm],
					   pba->a_today/a-1.0,
					   NULL,
					   &rho_ncdm,
					   &p_ncdm,
					   NULL,
					   NULL),
                   pba->error_message,
                   pba->error_message);
	rho_ncdm_rel_tot += 3.*p_ncdm;
	if (fabs(p_ncdm/rho_ncdm-1./3.)>ppr->tol_ncdm_initial_w)
	  is_early_enough = _FALSE_;
      }
      if (is_early_enough == _TRUE_)
	break;
      else
	a *= _SCALE_BACK_;
    }
    class_test(counter == _MAX_IT_,
	       pba->error_message,
	       "Search for initial scale factor a such that all ncdm species are relativistic failed.");
  }

  pvecback_integration[pba->index_bi_a] = a;

  /* Set initial values of {B} variables: */
  Omega_rad = pba->Omega0_g;
  if (pba->has_ur == _TRUE_)
    Omega_rad += pba->Omega0_ur;
  rho_rad = Omega_rad*pow(pba->H0,2)/pow(a/pba->a_today,4);
  if (pba->has_ncdm == _TRUE_){
    /** - We must add the relativistic contribution from NCDM species */
    rho_rad += rho_ncdm_rel_tot;
  }
  if (pba->has_dcdm == _TRUE_){
    /* Remember that the critical density today in CLASS conventions is H0^2 */
    pvecback_integration[pba->index_bi_rho_dcdm] =
      pba->Omega_ini_dcdm*pba->H0*pba->H0*pow(pba->a_today/a,3);
    if (pba->background_verbose > 3)
      printf("Density is %g. a_today=%g. Omega_ini=%g\n",pvecback_integration[pba->index_bi_rho_dcdm],pba->a_today,pba->Omega_ini_dcdm);
  }

  if (pba->has_dr == _TRUE_){
    if (pba->has_dcdm == _TRUE_){
      /**  - f is the critical density fraction of DR. The exact solution is:
       *
       * `f = -Omega_rad+pow(pow(Omega_rad,3./2.)+0.5*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3),2./3.);`
       *
       * but it is not numerically stable for very small f which is always the case.
       * Instead we use the Taylor expansion of this equation, which is equivalent to
       * ignoring f(a) in the Hubble rate.
       */
      f = 1./3.*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3)/sqrt(Omega_rad);
      pvecback_integration[pba->index_bi_rho_dr] = f*pba->H0*pba->H0/pow(a/pba->a_today,4);
    }
    else{
      /** There is also a space reserved for a future case where dr is not sourced by dcdm */
      pvecback_integration[pba->index_bi_rho_dr] = 0.0;
    }
  }

  if (pba->has_fld == _TRUE_){

    /* rho_fld today */
    rho_fld_today = pba->Omega0_fld * pow(pba->H0,2);

    /* integrate rho_fld(a) from a_ini to a_0, to get rho_fld(a_ini) given rho_fld(a0) */
    class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, pba->error_message);

    /* Note: for complicated w_fld(a) functions with no simple
    analytic integral, this is the place were you should compute
    numerically the simple 1d integral [int_{a_ini}^{a_0} 3
    [(1+w_fld)/a] da] (e.g. with the Romberg method?) instead of
    calling background_w_fld */

    /* rho_fld at initial time */
    pvecback_integration[pba->index_bi_rho_fld] = rho_fld_today * exp(integral_fld);

  }

  /** - Fix initial value of \f$ \phi, \phi' \f$
   * set directly in the radiation attractor => fixes the units in terms of rho_ur
   *
   * TODO:
   * - There seems to be some small oscillation when it starts.
   * - Check equations and signs. Sign of phi_prime?
   * - is rho_ur all there is early on?
   */
  if(pba->has_scf == _TRUE_){
    scf_lambda = pba->scf_parameters[0];
    if(pba->attractor_ic_scf == _TRUE_){
      pvecback_integration[pba->index_bi_phi_scf] = -1/scf_lambda*
        log(rho_rad*4./(3*pow(scf_lambda,2)-12))*pba->phi_ini_scf;
      if (3.*pow(scf_lambda,2)-12. < 0){
        /** - --> If there is no attractor solution for scf_lambda, assign some value. Otherwise would give a nan.*/
    	pvecback_integration[pba->index_bi_phi_scf] = 1./scf_lambda;//seems to the work
	if (pba->background_verbose > 0)
	  printf(" No attractor IC for lambda = %.3e ! \n ",scf_lambda);
      }
      pvecback_integration[pba->index_bi_phi_prime_scf] = 2*pvecback_integration[pba->index_bi_a]*
        sqrt(V_scf(pba,pvecback_integration[pba->index_bi_phi_scf]))*pba->phi_prime_ini_scf;
    }
    else{
      printf("Not using attractor initial conditions\n");
      /** - --> If no attractor initial conditions are assigned, gets the provided ones. */
      pvecback_integration[pba->index_bi_phi_scf] = pba->phi_ini_scf;
      pvecback_integration[pba->index_bi_phi_prime_scf] = pba->phi_prime_ini_scf;
    }
    class_test(!isfinite(pvecback_integration[pba->index_bi_phi_scf]) ||
               !isfinite(pvecback_integration[pba->index_bi_phi_scf]),
               pba->error_message,
               "initial phi = %e phi_prime = %e -> check initial conditions",
               pvecback_integration[pba->index_bi_phi_scf],
               pvecback_integration[pba->index_bi_phi_scf]);
  }

  /* Infer pvecback from pvecback_integration */
  class_call(background_functions(pba, pvecback_integration, pba->normal_info, pvecback),
	     pba->error_message,
	     pba->error_message);

  /* Just checking that our initial time indeed is deep enough in the radiation
     dominated regime */
  class_test(fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r,
	     pba->error_message,
	     "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
	     pvecback[pba->index_bg_Omega_r]);

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[pba->index_bg_H] <= 0.,
             pba->error_message,
             "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_tau] = 1./(a * pvecback[pba->index_bg_H]);

  /** - compute initial sound horizon, assuming \f$ c_s=1/\sqrt{3} \f$ initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau]/sqrt(3.);

  /** - set initial value of D and D' in RD. D will be renormalised later, but D' must be correct. */
  pvecback_integration[pba->index_bi_D] = a;
  pvecback_integration[pba->index_bi_D_prime] = 2*pvecback_integration[pba->index_bi_D]*pvecback[pba->index_bg_H];

  /** ---------------SELF INTERACTIONS BLOCK------------------*/

  if(pba->background_si_verbose>2) printf("%s\n", "Started executing SI block in background_initial_conditions()");

  if(pba->background_si_verbose>2){
    printf("Initial a = %e \n", a );
    printf("Initial H= %e \n", pvecback[pba->index_bg_H] );
  }

  double initial_reltime;

  //printf("Allocated initial_reltime\n");

  /* We first check if the species is coupled at the earliest computed time */

  if (pba->has_si_ncdm == _TRUE_){

    //printf("has_si_ncdm set to TRUE \n");

    for (int n_ncdm = 0; n_ncdm<pba->N_ncdm; n_ncdm++) {

      pba->is_early_coupled[n_ncdm] = _FALSE_;

      if(pba->ncdm_si_type[n_ncdm] >= 1){

        class_call(background_si_ncdm_reltime(pba,
          &initial_reltime,
          pba->a_today/a - 1.0,
          n_ncdm),
        pba->error_message,
        pba->error_message);

        if(pba->background_si_verbose>2)printf("Initial Relaxation time = %e \n", initial_reltime);

        if (initial_reltime <= 1./(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a])){
          if(pba->background_si_verbose>1){
            printf("Self Interacting NCDM species %d coupled at earliest computed time. Using equilibrium DF in background calculations.\n", n_ncdm+1);
          }
          pba->is_early_coupled[n_ncdm] = _TRUE_;
        }
        else {
          if(pba->background_si_verbose>0){
            printf("Self Interacting NCDM species %d is _NOT_ coupled at earliest computed time. Ignoring effects of SI.\n", n_ncdm + 1);
            printf("Remember! Effects of recoupling SI not implemented!\n");
          }
          pba->is_early_coupled[n_ncdm] = _FALSE_;
        }

        /* if early coupled, calculate equivalent parameters! */

        if (pba->is_early_coupled[n_ncdm] == _TRUE_){

          /** We had set a few flags in input.c to see if we need to infer particle mass or Omega.
           * Depending of how many of these were given, CLASS calculates m/Omega/normalization of f_0 to match DM abundance
           * We will need to do these things manually later if we re-initialize the SIWDM species*/
          if(pba->background_si_verbose>2){
            printf("mass_set_by_DF (before background_si_ncdm_equivalent_params) for species %d is set to %d\n", n_ncdm+1, pba->mass_set_by_DF[n_ncdm] );
            printf("Omega_set_by_DF (before background_si_ncdm_equivalent_params) for species %d is set to %d\n", n_ncdm+1, pba->Omega_set_by_DF[n_ncdm] );
          }

          /** So, the logic would be: if DF is set to FD then do nothing. If DF is passed by file, then:
          -> Calcualate its equivalent temperature
          -> If either m or Omega were infered, recalculate them
          -> If neither of them were infered, check that the factor (with the new temperature) is correct */

          if (pba->got_files[n_ncdm]){

            double equivalent_T;

            class_call(background_si_ncdm_equivalent_params(ppr, pba,
              &equivalent_T,
              NULL,
              pba->a_today/a-1.0,
              n_ncdm),
            pba->error_message,
            pba->error_message);

            if(pba->background_si_verbose>2){
              printf("mass_set_by_DF (after background_si_ncdm_equivalent_params) for species %d is set to %d\n", n_ncdm+1, pba->mass_set_by_DF[n_ncdm] );
              printf("Omega_set_by_DF (after background_si_ncdm_equivalent_params) for species %d is set to %d\n", n_ncdm+1, pba->Omega_set_by_DF[n_ncdm] );
            }

            pba->got_files[n_ncdm] = _FALSE_;
            pba->T_ncdm[n_ncdm] = equivalent_T;

            if(pba->background_si_verbose>2){
              printf("got_files for species %d is set to %d \n", n_ncdm + 1, pba->got_files[n_ncdm] );
              printf("T_ncdm for species %d is set to %.10e \n", n_ncdm + 1, pba->T_ncdm[n_ncdm] );
            }

          }

        }

      }

    }

    /* Here, after we have checked if the species has an initial coupled regime, we can check if the species remains self coupled
      into the Non Relativistic Regime...  */

    for (int n_ncdm=0; n_ncdm< pba->N_ncdm; n_ncdm++ ){

      pba->is_NR_SI_decoupled[n_ncdm] = _FALSE_;

      if (pba->ncdm_si_type[n_ncdm] >= 1){

        if(pba->background_si_verbose>2)printf("Checking if the ncdm species %d has coupled SI in the NR regime...\n", n_ncdm + 1);

        // a at T=m
        double a_nonrel = ( pba->T_cmb * pba->a_today * _k_B_ ) / (pba->m_ncdm_in_eV[n_ncdm] * _eV_);

        double * pvecback_integration_nonrel;
        double * pvecback_nonrel;
        class_alloc(pvecback_nonrel,pba->bg_size*sizeof(double),pba->error_message);
        class_alloc(pvecback_integration_nonrel,pba->bi_size*sizeof(double),pba->error_message);

        class_test(pba->has_fld || pba->has_dcdm || pba->has_scf || pba->has_dr, pba->error_message, 
          "Self Interacting NCDM currently incompatible with: UR, DCDM, SCF or DR.");

        // We quickly integrate the background until T=m
        pvecback_integration_nonrel[pba->index_bi_a] = a_nonrel;
        class_call(background_functions(pba, pvecback_integration_nonrel, pba->long_info, pvecback_nonrel),
          pba->error_message,
          pba->error_message);

        int n_si_ncdm = pba->ncdm_si_index[n_ncdm];

        printf("H_nonrel = %e, a_nonrel = %e, tau_nonrel=%e \n", pvecback_nonrel[pba->index_bg_H], pvecback_integration_nonrel[pba->index_bi_a], pvecback_nonrel[pba->index_bg_taurel_si_ncdm1+n_si_ncdm]);

        // Check if it is self-coupled
        if (pvecback_nonrel[pba->index_bg_taurel_si_ncdm1+n_si_ncdm]<=
          (1./(pvecback_nonrel[pba->index_bg_H]))){
          pba->is_NR_SI_decoupled[n_ncdm] = _TRUE_;
        }
        else{ 
          pba->is_NR_SI_decoupled[n_ncdm] = _FALSE_;
        }

        if(pba->background_si_verbose>2)printf("is_NR_SI_decoupled[n_ncdm] = %d\n", pba->is_NR_SI_decoupled[n_ncdm]);

      }

    }

    //printf("sum_early_coupled set to 0\n");

    if( pba->has_si_ncdm == _TRUE_ ){
      pba->sum_early_coupled=0;
      pba->sum_NR_SI_decoupled=0;
      for (n_ncdm=0;n_ncdm<pba->N_ncdm;n_ncdm++) {
        if(pba->is_early_coupled[n_ncdm] ==_TRUE_ && pba->ncdm_si_type >= 0 ) pba->sum_early_coupled++;
        if(pba->is_NR_SI_decoupled[n_ncdm] ==_TRUE_ && pba->ncdm_si_type >= 0) pba->sum_NR_SI_decoupled++;
      }
    }

    if(pba->background_si_verbose>0)printf("The total of early coupled species is %d\n",pba->sum_early_coupled);
    if(pba->background_si_verbose>0)printf("The total of NR self-decoupled species is %d\n",pba->sum_NR_SI_decoupled);

    /* On this block we perform the re-initialization of the ncdm parameters if it is early coupled and/or if it undergoes NR SI decoupling... */

    if((pba->sum_early_coupled>0)||(pba->sum_NR_SI_decoupled>0)){

      class_call(background_ncdm_init(ppr,pba),
       pba->error_message,
       pba->error_message); 

     /* Copied block from input.c -> Recalculate m, Omega, factor, etc.
      In this case, we will just use the flags set beforehand. */

     /* We must calculate M from omega or vice versa if one of them is missing.
       If both are present, we must update the degeneracy parameter to
       reflect the implicit normalization of the distribution function.*/

      double fnu_factor;
      double rho_ncdm;

      for (int n=0; n < pba->N_ncdm; n++){

        if (pba->ncdm_si_type[n] >= 0 && pba->is_early_coupled[n]== _TRUE_ && pba->got_files[n] == _TRUE_){

          printf("Species %d is early coupled and not set to thermal values. Recalculating factors for ncdm species %d \n", n+1, n+1 );

          pba->Omega0_ncdm_tot -= pba->Omega0_ncdm[n];
          if(pba->background_si_verbose>2)printf("Removing %g from total Omega..\n",pba->Omega0_ncdm[n]);

          if(pba->background_si_verbose>2){
            printf("The DF factor for species i=%d is set to %g.\n", n+1, pba->factor_ncdm[n] );
            printf("The mass value for species i=%d is set to %g eV.\n", n+1, pba->m_ncdm_in_eV[n] );
          }

          if (pba->mass_set_by_DF[n] == _FALSE_){
            /* Case of only mass or mass and Omega/omega: */
            pba->M_ncdm[n] = pba->m_ncdm_in_eV[n]/_k_B_*_eV_/pba->T_ncdm[n]/pba->T_cmb;
            class_call(background_ncdm_momenta(
             pba,
             n,
             pba->q_ncdm_bg[n],
             pba->w_ncdm_bg[n],
             pba->q_size_ncdm_bg[n],
             pba->M_ncdm[n],
             pba->factor_ncdm[n],
             0.,
             NULL,
             &rho_ncdm,
             NULL,
             NULL,
             NULL),
            pba->error_message,
            pba->error_message);
            if (pba->Omega_set_by_DF[n] == _TRUE_ ){
              pba->Omega0_ncdm[n] = rho_ncdm/pba->H0/pba->H0;
            }
            else if (pba->Omega_set_by_DF[n] == _FALSE_){
              fnu_factor = (pba->H0*pba->H0*pba->Omega0_ncdm[n]/rho_ncdm);
              pba->factor_ncdm[n] *= fnu_factor;
          /* dlnf0dlnq is already computed, but it is
             independent of any normalization of f0.
             We don't need the factor anymore, but we
             store it nevertheless:*/
              pba->deg_ncdm[n] *=fnu_factor;
            }
          }
          else if (pba->mass_set_by_DF[n] == _TRUE_){
            /* Case of only Omega/omega: */
            class_call(background_ncdm_M_from_Omega(ppr,pba,n),
             pba->error_message,
             pba->error_message);
            //printf("M_ncdm:%g\n",pba->M_ncdm[n]);
            pba->m_ncdm_in_eV[n] = _k_B_/_eV_*pba->T_ncdm[n]*pba->M_ncdm[n]*pba->T_cmb;
          }

          pba->Omega0_ncdm_tot += pba->Omega0_ncdm[n];
          if(pba->background_si_verbose>2)printf("Adding %g to total Omega..\n",pba->Omega0_ncdm[n]);

          if(pba->background_si_verbose>0){
            printf("The DF factor for species i=%d is set to %g.\n", n+1, pba->factor_ncdm[n] );
            printf("The mass value for species i=%d is set to %g eV.\n", n+1, pba->m_ncdm_in_eV[n] );
          }

        }

      }

    /** Provide the user with updated information about the ncdm species if some of these parameters are recalculated... */

      double Neff;
      double rho_ncdm_rel,rho_nu_rel;

      Neff = pba->Omega0_ur/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;

      /* loop over ncdm species */
      for (int n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

        /* call this function to get rho_ncdm */
        background_ncdm_momenta(
          pba,
          n_ncdm,
          pba->q_ncdm_bg[n_ncdm],
          pba->w_ncdm_bg[n_ncdm],
          pba->q_size_ncdm_bg[n_ncdm],
          0.,
          pba->factor_ncdm[n_ncdm],
          0.,
          NULL,
          &rho_ncdm_rel,
          NULL,
          NULL,
          NULL);

        /* inform user of the contribution of each species to
           radiation density (in relativistic limit): should be
           between 1.01 and 1.02 for each active neutrino species;
           evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
           density of one neutrino in the instantaneous decoupling
           limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
           from the definition of N_eff) */
        rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
        pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4);

        if(pba->background_si_verbose>0){
          printf("Recalculated sampling for ncdm species %d \n", n_ncdm );
          printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
           n_ncdm+1,
           pba->q_size_ncdm_bg[n_ncdm],
           pba->q_size_ncdm[n_ncdm],
           rho_ncdm_rel/rho_nu_rel);
        }

        Neff += rho_ncdm_rel/rho_nu_rel;

      }

      if(pba->background_si_verbose>0)printf(" -> total N_eff = %g (sumed over ultra-relativistic and ncdm species)\n",Neff);

    }

  }

  /** ------------END OF SELF INTERACTIONS BLOCK--------------*/

  return _SUCCESS_;

}

/**
 * Find the time of radiation/matter equality and store characteristic
 * quantitites at that time in the background structure..
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input/Output: pointer to background structure
  * @return the error status
 */

int background_find_equality(
                             struct precision *ppr,
                             struct background *pba) {

  double Omega_m_over_Omega_r=0.;
  int index_tau_minus = 0;
  int index_tau_plus = pba->bt_size-1;
  int index_tau_mid = 0;
  double tau_minus,tau_plus,tau_mid=0.;
  double * pvecback;

  /* first bracket the right tau value between two consecutive indices in the table */

  while ((index_tau_plus - index_tau_minus) > 1) {

    index_tau_mid = (int)(0.5*(index_tau_plus+index_tau_minus));

    Omega_m_over_Omega_r = pba->background_table[index_tau_mid*pba->bg_size+pba->index_bg_Omega_m]
      /pba->background_table[index_tau_mid*pba->bg_size+pba->index_bg_Omega_r];

    if (Omega_m_over_Omega_r > 1)
      index_tau_plus = index_tau_mid;
    else
      index_tau_minus = index_tau_mid;

  }

  /* then get a better estimate within this range */

  tau_minus = pba->tau_table[index_tau_minus];
  tau_plus =  pba->tau_table[index_tau_plus];

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  while ((tau_plus - tau_minus) > ppr->tol_tau_eq) {

    tau_mid = 0.5*(tau_plus+tau_minus);

    class_call(background_at_tau(pba,tau_mid,pba->long_info,pba->inter_closeby,&index_tau_minus,pvecback),
               pba->error_message,
               pba->error_message);

    Omega_m_over_Omega_r = pvecback[pba->index_bg_Omega_m]/pvecback[pba->index_bg_Omega_r];

    if (Omega_m_over_Omega_r > 1)
      tau_plus = tau_mid;
    else
      tau_minus = tau_mid;

  }

  pba->a_eq = pvecback[pba->index_bg_a];
  pba->H_eq = pvecback[pba->index_bg_H];
  pba->z_eq = pba->a_today/pba->a_eq -1.;
  pba->tau_eq = tau_mid;

  if (pba->background_verbose > 0) {
    printf(" -> radiation/matter equality at z = %f\n",pba->z_eq);
    printf("    corresponding to conformal time = %f Mpc\n",pba->tau_eq);
  }

  free(pvecback);

  return _SUCCESS_;

}


/**
 * Subroutine for formatting background output
 *
 */

int background_output_titles(struct background * pba,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){

  /** - Length of the column title should be less than _OUTPUTPRECISION_+6
      to be indented correctly, but it can be as long as . */
  int n;
  char tmp[24];
  int n_si_ncdm = 0;

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"comov. dist.",_TRUE_);
  class_store_columntitle(titles,"ang.diam.dist.",_TRUE_);
  class_store_columntitle(titles,"lum. dist.",_TRUE_);
  class_store_columntitle(titles,"comov.snd.hrz.",_TRUE_);
  class_store_columntitle(titles,"(.)rho_g",_TRUE_);
  class_store_columntitle(titles,"(.)rho_b",_TRUE_);
  class_store_columntitle(titles,"(.)rho_cdm",pba->has_cdm);
  if (pba->has_ncdm == _TRUE_){
    for (n=0; n<pba->N_ncdm; n++){
      sprintf(tmp,"(.)rho_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)p_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);

      /* SI */
      if (pba->ncdm_si_type[n]>0){
        if(pba->background_si_verbose>4)printf("Storing column titles...\n");
        n_si_ncdm++;
        sprintf(tmp,"(.)tau_ncdm[%d]",n);
        class_store_columntitle(titles,tmp,_TRUE_);
      }
      /* SI */
    }
  }
  class_store_columntitle(titles,"(.)rho_lambda",pba->has_lambda);
  class_store_columntitle(titles,"(.)rho_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)w_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)rho_ur",pba->has_ur);
  class_store_columntitle(titles,"(.)rho_crit",_TRUE_);
  class_store_columntitle(titles,"(.)rho_dcdm",pba->has_dcdm);
  class_store_columntitle(titles,"(.)rho_dr",pba->has_dr);

  class_store_columntitle(titles,"(.)rho_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_scf",pba->has_scf);
  class_store_columntitle(titles,"phi_scf",pba->has_scf);
  class_store_columntitle(titles,"phi'_scf",pba->has_scf);
  class_store_columntitle(titles,"V_scf",pba->has_scf);
  class_store_columntitle(titles,"V'_scf",pba->has_scf);
  class_store_columntitle(titles,"V''_scf",pba->has_scf);

  class_store_columntitle(titles,"gr.fac. D",_TRUE_);
  class_store_columntitle(titles,"gr.fac. f",_TRUE_);

  return _SUCCESS_;
}

int background_output_data(
                           struct background *pba,
                           int number_of_titles,
                           double *data){
  int index_tau, storeidx, n;
  double *dataptr, *pvecback;
  int n_si_ncdm;

  /** Stores quantities */
  for (index_tau=0; index_tau<pba->bt_size; index_tau++){
    dataptr = data + index_tau*number_of_titles;
    pvecback = pba->background_table + index_tau*pba->bg_size;
    storeidx = 0;

    class_store_double(dataptr,pba->a_today/pvecback[pba->index_bg_a]-1.,_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_time]/_Gyr_over_Mpc_,_TRUE_,storeidx);
    class_store_double(dataptr,pba->conformal_age-pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ang_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_lum_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rs],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_g],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_b],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_cdm],pba->has_cdm,storeidx);
    if (pba->has_ncdm == _TRUE_){
      n_si_ncdm = 0;
      for (n=0; n<pba->N_ncdm; n++){
        class_store_double(dataptr,pvecback[pba->index_bg_rho_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_p_ncdm1+n],_TRUE_,storeidx);

        /*SI*/
        if (pba->ncdm_si_type[n]>0){
          class_store_double(dataptr,pvecback[pba->index_bg_taurel_si_ncdm1+n_si_ncdm],_TRUE_,storeidx);
          n_si_ncdm++;
        } 
        /*SI*/
      }
    }
    class_store_double(dataptr,pvecback[pba->index_bg_rho_lambda],pba->has_lambda,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_w_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_ur],pba->has_ur,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_crit],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dcdm],pba->has_dcdm,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dr],pba->has_dr,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_rho_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_V_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_dV_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ddV_scf],pba->has_scf,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_D],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_f],_TRUE_,storeidx);
  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int background_derivs(
                      double tau,
                      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
                      double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  /** Summary: */

  /** - define local variables */

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback, a, H, rho_M;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - calculate functions of \f$ a \f$ with background_functions() */
  //printf("Called background_functions at background_derivs\n");
  class_call(background_functions(pba, y, pba->normal_info, pvecback),
             pba->error_message,
             error_message);

  /** - Short hand notation */
  a = y[pba->index_bi_a];
  H = pvecback[pba->index_bg_H];

  /** - calculate \f$ a'=a^2 H \f$ */
  dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

  /** - calculate \f$ t' = a \f$ */
  dy[pba->index_bi_time] = y[pba->index_bi_a];

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  /** - calculate \f$ rs' = c_s \f$*/
  dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]))*sqrt(1.-pba->K*y[pba->index_bi_rs]*y[pba->index_bi_rs]); // TBC: curvature correction

  /** - solve second order growth equation  \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */
  rho_M = pvecback[pba->index_bg_rho_b];
  if (pba->has_cdm)
    rho_M += pvecback[pba->index_bg_rho_cdm];
  dy[pba->index_bi_D] = y[pba->index_bi_D_prime];
  dy[pba->index_bi_D_prime] = -a*H*y[pba->index_bi_D_prime] + 1.5*a*a*rho_M*y[pba->index_bi_D];

  if (pba->has_dcdm == _TRUE_){
    /** - compute dcdm density \f$ \rho' = -3aH \rho - a \Gamma \rho \f$*/
    dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dcdm]-
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if ((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_)){
    /** - Compute dr density \f$ \rho' = -4aH \rho - a \Gamma \rho \f$ */
    dy[pba->index_bi_rho_dr] = -4.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dr]+
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if (pba->has_fld == _TRUE_) {
    /** - Compute fld density \f$ \rho' = -3aH (1+w_{fld}(a)) \rho \f$ */
    dy[pba->index_bi_rho_fld] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*(1.+pvecback[pba->index_bg_w_fld])*y[pba->index_bi_rho_fld];
  }

  if (pba->has_scf == _TRUE_){
    /** - Scalar field equation: \f$ \phi'' + 2 a H \phi' + a^2 dV = 0 \f$  (note H is wrt cosmic time) */
    dy[pba->index_bi_phi_scf] = y[pba->index_bi_phi_prime_scf];
    dy[pba->index_bi_phi_prime_scf] = - y[pba->index_bi_a]*
      (2*pvecback[pba->index_bg_H]*y[pba->index_bi_phi_prime_scf]
       + y[pba->index_bi_a]*dV_scf(pba,y[pba->index_bi_phi_scf])) ;
  }


  return _SUCCESS_;

}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 * For Albrecht & Skordis model: 9908085
 * - \f$ V = V_{p_{scf}}*V_{e_{scf}} \f$
 * - \f$ V_e =  \exp(-\lambda \phi) \f$ (exponential)
 * - \f$ V_p = (\phi - B)^\alpha + A \f$ (polynomial bump)
 *
 * TODO:
 * - Add some functionality to include different models/potentials (tuning would be difficult, though)
 * - Generalize to Kessence/Horndeski/PPF and/or couplings
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * - Numerical derivatives may further serve as a consistency check.
 *
 */

/**
 *
 * The units of phi, tau in the derivatives and the potential V are the following:
 * - phi is given in units of the reduced Planck mass \f$ m_{pl} = (8 \pi G)^{(-1/2)}\f$
 * - tau in the derivative is given in units of Mpc.
 * - the potential \f$ V(\phi) \f$ is given in units of \f$ m_{pl}^2/Mpc^2 \f$.
 * With this convention, we have
 * \f$ \rho^{class} = (8 \pi G)/3 \rho^{physical} = 1/(3 m_{pl}^2) \rho^{physical} = 1/3 * [ 1/(2a^2) (\phi')^2 + V(\phi) ] \f$
    and \f$ \rho^{class} \f$ has the proper dimension \f$ Mpc^-2 \f$.
 */

double V_e_scf(struct background *pba,
               double phi
               ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return  exp(-scf_lambda*phi);
}

double dV_e_scf(struct background *pba,
                double phi
                ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return -scf_lambda*V_scf(pba,phi);
}

double ddV_e_scf(struct background *pba,
                 double phi
                 ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return pow(-scf_lambda,2)*V_scf(pba,phi);
}


/** parameters and functions for the polynomial coefficient
 * \f$ V_p = (\phi - B)^\alpha + A \f$(polynomial bump)
 *
 * double scf_alpha = 2;
 *
 * double scf_B = 34.8;
 *
 * double scf_A = 0.01; (values for their Figure 2)
 */

double V_p_scf(
               struct background *pba,
               double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double dV_p_scf(
                struct background *pba,
                double phi) {

  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return   scf_alpha*pow(phi -  scf_B,  scf_alpha - 1);
}

double ddV_p_scf(
                 struct background *pba,
                 double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  scf_alpha*(scf_alpha - 1.)*pow(phi -  scf_B,  scf_alpha - 2);
}

/** Fianlly we can obtain the overall potential \f$ V = V_p*V_e \f$
 */

double V_scf(
             struct background *pba,
             double phi) {
  return  V_e_scf(pba,phi)*V_p_scf(pba,phi);
}

double dV_scf(
              struct background *pba,
	      double phi) {
  return dV_e_scf(pba,phi)*V_p_scf(pba,phi) + V_e_scf(pba,phi)*dV_p_scf(pba,phi);
}

double ddV_scf(
               struct background *pba,
               double phi) {
  return ddV_e_scf(pba,phi)*V_p_scf(pba,phi) + 2*dV_e_scf(pba,phi)*dV_p_scf(pba,phi) + V_e_scf(pba,phi)*ddV_p_scf(pba,phi);
}
