#A simple python script to compare SIWDM, WDM and LCDM power spectra
#Based on the original implementation of scripts/warmup.py in CLASS 2.7
# coding: utf-8

# In[ ]:

# import classy module
from classy import Class

# In[ ]:

# define functions to get Pk, cls outputs from class "Class"
def Class_getcls(Class):
    # get all C_l output
    cls = Class.lensed_cl(2500)
    # To check the format of cls
    cls.viewkeys()

    ll = cls['ell'][2:]
    clTT = cls['tt'][2:]
    clEE = cls['ee'][2:]
    clPP = cls['pp'][2:]

    return {'ll': ll, 'clTT': clTT, 'clEE': clEE, 'clPP': clPP}

import numpy as np
def Class_getpk(Class):
    # get P(k) at redhsift z=0
    # The function .pk(k,z) wants k in 1/Mpc so we must convert k for each case with the right h
    h = LambdaCDM.h() # get reduced Hubble for conversions to 1/Mpc
    return lambda k_vec: np.array([Class.pk(k*h,0.)*h**3 for k in k_vec]) #function .pk returns in Mpc^3 so we must convert to [Mpc/h]^3


# In[ ]:


# define common parameters
common_settings = {'omega_b':0.022032,
'h':0.67556,
'A_s':2.215e-9,
'n_s':0.9619,
'tau_reio':0.0925, 
'output':'tCl,pCl,lCl,mPk,mTk', 
'lensing':'yes', 
'P_k_max_1/Mpc':300., 
'perturbations_verbose': 3}

# create instance of the class "Class" for LCDM
LambdaCDM = Class()
#pass LambdaCDM parameters
LambdaCDM.set(common_settings)
LambdaCDM.set({'omega_cdm':0.12038})
# run class
LambdaCDM.compute()

# # In[ ]:

# create instance of the class "Class" for WDM
WDM = Class()
#pass LambdaCDM parameters
WDM.set(common_settings)
#set input parameters
WDM.set({'omega_cdm': 0.0000001, 
    'N_ncdm': 1,
    'use_ncdm_psd_files': 0,
    'm_ncdm': 1.111e+04,
    'omega_ncdm': 0.1201074,
    'T_ncdm': 7.13766e-01,
    'ncdm_self_interaction_type': 0})
#set precision parameters
WDM.set({# precision file for a quick run of expanatory_plus_siwdm.ini
        'hyper_flat_approximation_nu' : 7000.,
        'transfer_neglect_delta_k_S_t0' : 0.17,
        'transfer_neglect_delta_k_S_t1' : 0.05,
        'transfer_neglect_delta_k_S_t2' : 0.17,
        'transfer_neglect_delta_k_S_e' : 0.13,
        'delta_l_max' : 1000,
        #Non stiff evolver
        'evolver' : 0,
        #Change if neccesary!
        'k_per_decade_for_pk' : 20,
        #Higher numbers (up to 40ish) can improve precision at higher k
        'l_max_ncdm' : 6,
        #Ignore ncdm fluid approx., important!
        'ncdm_fluid_approximation' : 3,
        #This is one of the most important precision parameters in this case, up to 1e-14 can improve precision
        'tol_ncdm_bg' : 1.000e-03,
        'tol_ncdm_synchronous' : 1.000e-02,
        'tol_ncdm_newtonian' : 1.000e-02,
        #Added params for siwdm! 
        #Estimation method for l>2 moments in tca approx. (0: \Psi_l>2 = 0, 1: tridiagonal matrix method) 
        'ncdm_ic_estimation_method' : 0,
        #Closing approxmination for hierarchies (0: Ma&Bertschiner method, 1:Unused)
        'ncdm_lmax_method' : 0,
        #Trigger for TCA approximation (in units of tau_rel/(tau_k or tau_h)) 
        'ncdm_si_tca_trigger' : 1.0})
# run class
WDM.compute()

# In[ ]:
#tell class where the SIWDM relaxation time table is defined
import os
reltime_path = os.getcwd() + '/explanatory_plus_siwdm_reltime.dat'

# create instance of the class "Class" for SIWDM
SIWDM = Class()
#pass LambdaCDM parameters
SIWDM.set(common_settings)
#set input parameters
SIWDM.set({'omega_cdm': 0.0000001, 
    'N_ncdm': 1,
    'use_ncdm_psd_files': 0,
    'm_ncdm': 1.111e+04,
    'omega_ncdm': 0.1201074,
    'T_ncdm': 7.13766e-01,
    'ncdm_self_interaction_type' : 1,
    'ncdm_self_interaction_tau_filenames' : reltime_path,
    'ncdm_NR_SI_decoupling_method' : 1,
    'background_si_verbose' : 1})
#set precision parameters
SIWDM.set({# precision file for a quick run of expanatory_plus_siwdm.ini
        'hyper_flat_approximation_nu' : 7000.,
        'transfer_neglect_delta_k_S_t0' : 0.17,
        'transfer_neglect_delta_k_S_t1' : 0.05,
        'transfer_neglect_delta_k_S_t2' : 0.17,
        'transfer_neglect_delta_k_S_e' : 0.13,
        'delta_l_max' : 1000,
        #Non stiff evolver
        'evolver' : 0,
        #Change if neccesary!
        'k_per_decade_for_pk' : 20,
        #Higher numbers (up to 40ish) can improve precision at higher k
        'l_max_ncdm' : 6,
        #Ignore ncdm fluid approx., important!
        'ncdm_fluid_approximation' : 3,
        #This is one of the most important precision parameters in this case, up to 1e-14 can improve precision
        'tol_ncdm_bg' : 1.000e-03,
        'tol_ncdm_synchronous' : 1.000e-02,
        'tol_ncdm_newtonian' : 1.000e-02,
        #Added params for siwdm! 
        #Estimation method for l>2 moments in tca approx. (0: \Psi_l>2 = 0, 1: tridiagonal matrix method) 
        'ncdm_ic_estimation_method' : 0,
        #Closing approxmination for hierarchies (0: Ma&Bertschiner method, 1:Unused)
        'ncdm_lmax_method' : 0,
        #Trigger for TCA approximation (in units of tau_rel/(tau_k or tau_h)) 
        'ncdm_si_tca_trigger' : 1.0})
# run class
SIWDM.compute()

# In[ ]:

# uncomment to get plots displayed in notebook
#%matplotlib inline
import matplotlib.pyplot as plt
from math import pi


# In[ ]:

# plot C_l^TT
plt.figure(1)
plt.xscale('log');plt.yscale('linear');plt.xlim(2,2500)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{TT}$')

#LambdaCDM
ll = Class_getcls(LambdaCDM)['ll']
clTT = Class_getcls(LambdaCDM)['clTT']
plt.plot(ll,clTT*ll*(ll+1)/2./pi,'r-', label= 'LambdaCDM')
#WDM
ll = Class_getcls(WDM)['ll']
clTT = Class_getcls(WDM)['clTT']
plt.plot(ll,clTT*ll*(ll+1)/2./pi,'b-', label= 'WDM')

#WDM
ll = Class_getcls(SIWDM)['ll']
clTT = Class_getcls(SIWDM)['clTT']
plt.plot(ll,clTT*ll*(ll+1)/2./pi,'g-', label= 'SIWDM')

plt.legend()


# In[ ]:

plt.savefig('Simple_SIWDM_cltt.pdf')


# In[ ]:

LambdaCDM_pk = Class_getpk(LambdaCDM)
WDM_pk = Class_getpk(WDM)
SIWDM_pk = Class_getpk(SIWDM)

# In[ ]:

# plot P(k)
kk = np.logspace(-4, np.log10(300.), 1000)
plt.figure(2)
plt.xscale('log');plt.yscale('log');plt.xlim(kk[0],kk[-1])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
plt.plot(kk,LambdaCDM_pk(kk),'r-', label= 'LambdaCDM')
plt.plot(kk,WDM_pk(kk),'g-', label= 'WDM')
plt.plot(kk,SIWDM_pk(kk),'b-', label= 'SIWDM')

plt.legend()


# In[ ]:

plt.savefig('Simple_SIWDM_pk.pdf')


# In[ ]:

# optional: clear content of LambdaCDM (to reuse it for another model)
LambdaCDM.struct_cleanup()
WDM.struct_cleanup()
SIWDM.struct_cleanup()
# optional: reset parameters to default
LambdaCDM.empty()
WDM.empty()
SIWDM.empty()

# In[ ]:
