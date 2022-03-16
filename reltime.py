'''
This simple script calculates relaxation time tables for SIWDM components. 
When using the SIWDM module, you need to supply the integrator with a time-ordered table with the values for the relaxation time
(as defined in Yunis et. al. 2022). In this case, what needs to be supplied is a table with the CMB temperature in eV versus
the relaxation time in Mpc. This script does just that, when it is supplied with the following information:

arguments:
- Particle Mass (mandatory): the SIWDM particle's mass, in eV
- Self Interaction Model (mandatory): the tree level process responsible for the collisions. Currently supported ones are (see Yunis et. al. 2020):
    Massive Mediator (Cm)
    Massive Vector Field (CV)
    Constant Amplitude (C0)
- Interaction Constant (mandatory): the coupling constant of the interaction model chosen, proportional to the cross section.
    For the Vector Field and Massive Mediator models, it's defined roughly as ( g^4 / (2 m_phi^4) in eV^(-4) )
    (for reference, for these models 1e-50 is relativistic decoupling, and 1e-30 is roughly the BC upper limit), 
    for the constant amplitude one it's 6g^4,
    with g the couplign constant
- Background model (mandatory): Roughly which combination of DM temperature and background distribution to use. Supported ones are:
    Sterile Neutrino Non-Resonant Production (nonresonant): T_DM ~ T_neutrino, f_DM ~ Fermi-Dirac, rescaled to the right DM abundance 
    Sterine Neutrino Resonant Production (steriledm): T_DM ~ T_neutrino, f_DM is externally supplied (check Venumadhav et. al. 2015 for an application to get these)
    Thermal Priduction (thermal): T_DM << T_neutrino, f_DM ~ Fermi-Dirac (not rescaled, just produced very early)
- Distribution Function File (only if steriledm is selected): path to the distribution function required to run resonant production (a table of q/T_ncdm vs f).

Typical call of this script is done by

python reltime.py [Particle Mass] --[Self Interaction Model] [Self Interaction Constant] --[Background Model] [Distribution Function File]

E.g.
python reltime.py 1.1e4 --CV 1e-30 --nonresonant

Would return a table roughly like the one in explanatory_plus_siwdm_reltime.dat

For the output of this script:
It will generate the table in %working_dir%/reltime/relaxation_time_table_n.dat, with a table on descending order in T_gamma. The script itself outputs
both its error code (first) and the file name where it wrote the data (second)
'''

import numpy as np
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy import interpolate
import matplotlib.pyplot as plt
import warnings

import sys
import math
import os
import shutil
import subprocess
import glob

def return_function(success_code,reltime_file):
    print(success_code,',',reltime_file)

def error_handler(error_code):
    sys.stderr.write('Exception ocurred in interface_reltime\n')
    if (error_code==-1):
        sys.stderr.write('Error code: -1, insufficient parameters. See documentation\n')
        return_function(-1,'')
        return -1
    elif (error_code==-2):
        sys.stderr.write('Error code: -2, Non regular numbers for mass or amplitude constant\n')
        return_function(-2,'')
        return -1
    elif (error_code==-3):
        sys.stderr.write('Error code: -3, Unidentified model type. Use one of --CV, --Cm, --C0\n')
        return_function(-3,'')
        return -1
    elif (error_code==-4):
        sys.stderr.write('Error code: -4, Unidentified psd type. Use one of --nonresonant, --steriledm, --thermal\n')
        return_function(-4,'')
        return -1
    elif (error_code==-5):
        sys.stderr.write('Error code: -5, distribution function file not provided for mode --steriledm\n')
        return_function(-5,'')
        return -1

def get_Tncdm0_reldec(model,psd_file,mass):
    if model == '--nonresonant':
        return np.power(4/11,1/3) * 0.000234822252 # approximately T_nu(active) (in eV)
    elif model == '--thermal':
        m = mass * 1e6 #mass in eV
        Omega_ncdm = 0.1201074
        Tx_over_Tnu = np.power(Omega_ncdm/m*94,1/3) # formula Viel 2005 (astro-ph/0501562)
        return Tx_over_Tnu * np.power(4/11,1/3) * 0.000234822252 # approximately T_nu(active) (in eV)
    elif model == '--steriledm':
        #Import data
        print(psd_file)
        data = np.loadtxt(psd_file)
        q = data[:,0]
        f = data[:,1]

        #define interpolator for f0
        tck = interpolate.splrep(q, f, s=0)
        def f_0(x):
            if isinstance(x, float):
                if x<min(q):
                    return f[0]
                if x>max(q):
                    return 0
                else:
                    return interpolate.splev(x, tck, der=0)
                
            if isinstance(x, np.ndarray):
                temp_f = np.zeros(x.size)

                temp_f[x<min(q)] = f[0]
                temp_f[x>max(q)] = 0
                temp_f[(x<max(q)) & (x>min(q))] = interpolate.splev(x[(x<max(q)) & (x>min(q))], tck, der=0)

                return temp_f

        #Integrate n, rho
        def integrand_n(q):
            return q**2 * f_0(q)

        #In this case, epsilon = sqrt(q^2+m^2a^2/T_ncdm^2): In this scenario this last term is m/T_ncdm(10MeV) = m/10MeV
        def integrand_rho(q):
            return q**2 * f_0(q) * np.sqrt( q**2 + mass**2/10**2 )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            n = 4 * np.pi * np.power(10, 3) * quad(integrand_n, 0, np.inf)[0]
            rho = 4 * np.pi * np.power(10, 4) * quad(integrand_rho, 0, np.inf)[0]

        # F_2/3(z=0) the complete fermi dirac integral of order 2/3 with 0 fugacity
        K_3 = 5.6822
        K_2 = 1.8030

        #calculate parameters
        k_bT = rho / n * K_2 / K_3
        C = rho / (k_bT**4) / K_3

        #return
        return (k_bT/10) * np.power(4/11,1/3) * 0.000234822252 # approximately (T_sterile/T_nu(active))|_10MeV * T_nu(active)|_today (in eV)

def get_normalization_factors(T_ncdm0, mass):

    m = mass * 1e6 #mass in eV

    rho_0_ncdm_reldec = (3/4) * 1.202 / (np.pi**2) * 2 * T_ncdm0**3 * m
    rho_crit = (2.77536627e11 #h^2 M_sol Mpc^-3
            * 0.7**2 #h^2
            * 1.11574e66 #M_sol in eV (Google)
            * np.power(1.5637e29,-3) ) #Mpc in eV^(-1) (Google)
    
    norm_reldec = 0.315*rho_crit/rho_0_ncdm_reldec

    T_ncdmNR0= 0.4652 * (T_ncdm0**2/m) #As from before. It is the same as if m * (aT/a0)**2 apart from the normalization factor!

    rho_0_ncdm_nrdec = ( m *
                    2 * #g
                    np.power(m*T_ncdmNR0/(2*np.pi), 3/2) )
                    
    norm_nrdec = 0.315*rho_crit/rho_0_ncdm_nrdec

    return norm_reldec, norm_nrdec

T_cmb0 = 0.000234822252 # approximately T_CBM (in eV)

def T_ncdm_reldec(a, m, T_ncdm0):
    return T_ncdm0 / a

def T_ncdm_nrdec(a, m, T_ncdm0):
    
    aT = T_ncdm0/m
    
    if type(a) is float:
        if (a<aT):
            return T_ncdm_reldec(a, m, T_ncdm0)
        if (a>aT):
            return 0.4652 * (T_ncdm0**2/m) / a**2
    if type(a) is np.ndarray:
        tempArray = np.zeros_like(a)
        tempArray[a<aT] = T_ncdm_reldec(a[a<aT], m, T_ncdm0)
        tempArray[a>aT] = 0.4652 * (T_ncdm0**2/m) * np.power(a[a>aT], -2)
        return tempArray

def n_ncdm_nr_nrdec(a, m, nrdec_factor, T_ncdm0):
    return 2 * nrdec_factor * np.power(m*T_ncdm_nrdec(a, m, T_ncdm0)/2/np.pi, 3/2)
def n_ncdm_rel_nrdec(a, m, nrdec_factor, T_ncdm0):
    return (7/8)*(1.202/np.pi**2) * 2 * (nrdec_factor/4.534) * np.power(T_ncdm_nrdec(a, m, T_ncdm0),3)
def n_ncdm_reldec(a, m, reldec_factor, T_ncdm0):
    return (7/8)*(1.202/np.pi**2) * 2 * reldec_factor * np.power(T_ncdm_reldec(a, m, T_ncdm0),3)

from scipy.special import kn

def GG_Gamma_mC0_norm(chi,T):
    
    #El input T es \tilde{T} \equiv T/m
    #Primero, la normalizacion de esta integral en s.
    
    Norm = ( 4 * np.pi * T * kn(2, 1/T) )
    
    #La integral en s
    integrand = lambda s: chi(s)*np.sqrt(1-4/s)*(np.sqrt(s)*T)*kn(1,np.sqrt(s)/T)
    
    #En general cae rapidamente. Metodo de integracion logaritmico!
    s = lambda b: np.exp(b)
    logIntegrand = lambda b: np.log(integrand(s(b)))
    fullIntegrand = lambda b: np.exp( logIntegrand(b) + b )
    
    #Valor maximo de la integral...
    Int_upper_limit = max( (20*T)**2, 5 )
    
    #Hago la integral
    Int_Tau_Val = quad(lambda b: fullIntegrand(b) , 
                       np.log(4) + 10*np.finfo(float).eps , 
                       np.log(Int_upper_limit) ) #+ 20*np.finfo(float).eps)
    
    #Devuelvo el valor de la integral
    Int_Tau = Int_Tau_Val[0]
    if (Norm != 0.0):
        invTauThAvg = Int_Tau/Norm
    else:
        invTauThAvg = float('inf')
        sys.stderr.write("Integral is thermally supressed")
    
    Gamma_mCo = invTauThAvg  #Esto me da \Gamma / C_i*m, sin unidades
    
    return Gamma_mCo

def Gamma_0_nrdec_NR(a, m, C0, T_ncdm0, nrdec_factor):
    T = T_ncdm_nrdec(a, m, T_ncdm0)
    sigma_v_0_nrdec= C0 / (4*np.sqrt(np.pi)*m**2) * ( np.power(T/m, 1/2) - np.power(T/m,3/2) + 9/4*np.power(T/m,5/2) )
    return sigma_v_0_nrdec * n_ncdm_nr_nrdec(a, m, nrdec_factor, T_ncdm0)

def Gamma_0_rdec_NR(a, m, C0, T_ncdm0, reldec_factor):
    T = T_ncdm_reldec(a, m , T_ncdm0)
    sigma_v_0_rdec = C0 * 35 / (1024*m**2) * ( 16*np.power(T/m, 1) - 72 * np.power(T/m,3) + 891 * np.power(T/m,5) )
    return sigma_v_0_rdec * n_ncdm_reldec(a, m, reldec_factor, T_ncdm0)

def Gamma_0_nrdec_R(a, m, C0, T_ncdm0, nrdec_factor):
    result = np.zeros_like(a)
    for i in range(0,a.size):
        result[i] = GG_Gamma_mC0_norm(chi_0, T_ncdm_reldec(a[i], m, T_ncdm0)/m) * m * C0 * (nrdec_factor/4.534)
    return result

def Gamma_0_rdec_R(a, m, C0, T_ncdm0, reldec_factor):
    result = np.zeros_like(a)
    for i in range(0,a.size):
        result[i] = GG_Gamma_mC0_norm(chi_0, T_ncdm_reldec(a[i], m, T_ncdm0)/m) * m * C0 * reldec_factor
    return result

def Gamma_m_nrdec_NR(a, m, Cm, T_ncdm0, nrdec_factor): #Remember to use Cm_over_m4
    T = T_ncdm_nrdec(a, m, T_ncdm0)
    sigma_v_m_nrdec= 4 * Cm * m**2 / np.sqrt(np.pi) * ( np.power(T/m, 1/2) + 3*np.power(T/m,3/2) + 137/4*np.power(T/m,5/2) )
    return sigma_v_m_nrdec * n_ncdm_nr_nrdec(a, m, nrdec_factor, T_ncdm0)

def Gamma_m_rdec_NR(a, m, Cm, T_ncdm0, reldec_factor): #Remember to use Cm_over_m4
    T = T_ncdm_reldec(a, m, T_ncdm0)
    sigma_v_m_rdec = 35/64 * Cm * m**2 * ( 16*np.power(T/m, 1) + 216 * np.power(T/m,3) + 13563 * np.power(T/m,5) )
    return sigma_v_m_rdec * n_ncdm_reldec(a, m, reldec_factor, T_ncdm0)

def Gamma_m_nrdec_R(a, m, Cm, T_ncdm0, nrdec_factor):
    result = np.zeros_like(a)
    for i in range(0,a.size):
        result[i] = GG_Gamma_mC0_norm(chi_m, T_ncdm_reldec(a[i], m, T_ncdm0)/m) * m * Cm * (nrdec_factor/4.534)
    return result

def Gamma_m_rdec_R(a, m, Cm, T_ncdm0, reldec_factor):
    result = np.zeros_like(a)
    for i in range(0,a.size):
        result[i] = GG_Gamma_mC0_norm(chi_m, T_ncdm_reldec(a[i], m, T_ncdm0)/m) * m * Cm * reldec_factor
    return result

def Gamma_V_nrdec_NR(a, m, CV, T_ncdm0, nrdec_factor):
    T = T_ncdm_nrdec(a, m, T_ncdm0)
    sigma_v_V_nrdec= CV * m**2 / np.sqrt(np.pi) * ( 42*np.power(T/m, 1/2) + 302*np.power(T/m,3/2) + 1013/2*np.power(T/m,5/2) )
    return sigma_v_V_nrdec * n_ncdm_nr_nrdec(a, m, nrdec_factor, T_ncdm0)

def Gamma_V_rdec_NR(a, m, CV, T_ncdm0, reldec_factor):
    T = T_ncdm_reldec(a, m, T_ncdm0)
    sigma_v_V_rdec = 105/128 * CV * m**2 * ( 112*np.power(T/m, 1) + 3624 * np.power(T/m,3) + 33429 * np.power(T/m,5) )
    return sigma_v_V_rdec * n_ncdm_reldec(a, m, reldec_factor, T_ncdm0)

def Gamma_V_nrdec_R(a, m, CV, T_ncdm0, nrdec_factor):
    result = np.zeros_like(a)
    for i in range(0,a.size):
        result[i] = GG_Gamma_mC0_norm(chi_V, T_ncdm_reldec(a[i], m, T_ncdm0)/m) * m * CV * (nrdec_factor/4.534)
    return result

def Gamma_V_rdec_R(a, m, CV, T_ncdm0, reldec_factor):
    result = np.zeros_like(a)
    for i in range(0,a.size):
        result[i] = GG_Gamma_mC0_norm(chi_V, T_ncdm_reldec(a[i], m, T_ncdm0)/m) * m * CV * reldec_factor
    return result

def get_R_slope(Gamma_R, m, C, T_ncdm0, reldec_factor):

    aT = T_ncdm0/m

    a_rel_1 = 1e-5 * aT
    a_rel_2 = 1e-3 * aT
    
    temp = Gamma_R(np.array([a_rel_1, a_rel_2]), m, C, T_ncdm0, reldec_factor)
    Gamma_rel_1=temp[0]
    Gamma_rel_2=temp[1]
    
    beta = (np.log(Gamma_rel_1) - np.log(Gamma_rel_2)) / (np.log(a_rel_1)-np.log(a_rel_2))
    alpha = np.exp( np.log(Gamma_rel_1) - beta * np.log(a_rel_1) )
    
    return np.array([alpha,beta])

def get_NR_slope(Gamma_NR, m, C, T_ncdm0, nrdec_factor):

    aT = T_ncdm0/m

    a_nrel_1 = 1e3 * aT
    a_nrel_2 = 1e5 * aT
    
    temp = Gamma_NR(np.array([a_nrel_1, a_nrel_2]), m, C, T_ncdm0, nrdec_factor)
    Gamma_nrel_1=temp[0]
    Gamma_nrel_2=temp[1]
    
    beta = (np.log(Gamma_nrel_1) - np.log(Gamma_nrel_2)) / (np.log(a_nrel_1)-np.log(a_nrel_2))
    alpha = np.exp( np.log(Gamma_nrel_1) - beta * np.log(a_nrel_1) )
    
    return np.array([alpha,beta])

def Gamma_0_nrdec_interp(a, m, C0, T_ncdm0, nrdec_factor):
    params_R = get_R_slope(Gamma_0_nrdec_R, m, C0, T_ncdm0, nrdec_factor)
    params_NR = get_NR_slope(Gamma_0_nrdec_NR, m, C0, T_ncdm0, nrdec_factor)
    return np.power( 0.5 * (np.power(params_NR[0]*np.power(a,params_NR[1]),-1) 
                            + np.power(params_R[0]*np.power(a,params_R[1]), -1)), -1)

def Gamma_0_rdec_interp(a, m, C0, T_ncdm0, reldec_factor):
    params_R = get_R_slope(Gamma_0_rdec_R, m, C0, T_ncdm0, reldec_factor)
    params_NR = get_NR_slope(Gamma_0_rdec_NR, m, C0, T_ncdm0, reldec_factor)
    return np.power( 0.5 * (np.power(params_NR[0]*np.power(a,params_NR[1]),-1)
                            + np.power(params_R[0]*np.power(a,params_R[1]), -1)), -1)

def Gamma_m_nrdec_interp(a, m, Cm, T_ncdm0, nrdec_factor):
    Cm_over_m4 = Cm/m**4
    params_R = get_R_slope(Gamma_m_nrdec_R, m, Cm, T_ncdm0, nrdec_factor)
    params_NR = get_NR_slope(Gamma_m_nrdec_NR, m, Cm_over_m4, T_ncdm0, nrdec_factor)
    return 0.5 * ( params_R[0]*np.power(a,params_R[1]) + params_NR[0]*np.power(a,params_NR[1]) )

def Gamma_m_rdec_interp(a, m, Cm, T_ncdm0, reldec_factor):
    Cm_over_m4 = Cm/m**4
    params_R = get_R_slope(Gamma_m_rdec_R, m, Cm, T_ncdm0, reldec_factor)
    params_NR = get_NR_slope(Gamma_m_rdec_NR, m, Cm_over_m4, T_ncdm0, reldec_factor)
    return 0.5 * ( params_R[0]*np.power(a,params_R[1]) + params_NR[0]*np.power(a,params_NR[1]) )

def Gamma_V_nrdec_interp(a, m, CV, T_ncdm0, nrdec_factor):
    CV_over_m4 = CV/m**4
    params_R = get_R_slope(Gamma_V_nrdec_R, m, CV, T_ncdm0, nrdec_factor)
    params_NR = get_NR_slope(Gamma_V_nrdec_NR, m, CV_over_m4, T_ncdm0, nrdec_factor)
    return 0.5 * ( params_R[0]*np.power(a,params_R[1]) + params_NR[0]*np.power(a,params_NR[1]) )

def Gamma_V_rdec_interp(a, m, CV, T_ncdm0, reldec_factor):
    CV_over_m4 = CV/m**4
    params_R = get_R_slope(Gamma_V_rdec_R, m, CV, T_ncdm0, reldec_factor)
    params_NR = get_NR_slope(Gamma_V_rdec_NR, m, CV_over_m4, T_ncdm0, reldec_factor)
    return 0.5 * ( params_R[0]*np.power(a,params_R[1]) + params_NR[0]*np.power(a,params_NR[1]) )


def chi_0(s):
    return np.sqrt(1-4/s)

def chi_m(s):
    return np.sqrt(1-4/s) * (256 - 128*s + 19*(s**2))

def chi_V(s):
    return np.sqrt(1-4/s) * (74 - 29*s) * (1 - s)

def unitConversion_eV_2_K(T_eV):
    return 1.16e4 * T_eV
def unitConversion_eVm1_2_seconds(t_eVm1):
    return 6.58e-16 * t_eVm1
def unitConversion_time_evm1_2_Mpc(t_eVm1):
    longitude_in_meters = unitConversion_eVm1_2_seconds(t_eVm1) * 299792458. # time*c
    return longitude_in_meters * 3.24e-23 #Megaparsecs

def H(a):
# Data from http://pdg.lbl.gov/2019/reviews/rpp2019-rev-cosmological-parameters.pdf and Planck 2018
    return 5.6759e-32 * np.sqrt( 5.502e-5 * np.power(a,-4) + 0.315 * np.power(a,-3) ) # HO*sqrt(Omega_0r*a^-4+Omega_0m*a^-3)

def is_reldec(Gamma_rdec_interp, Gamma_nrdec_interp, m, C, T_ncdm0, reldec_factor):
    a_equality = 1./(3601.)
    if (Gamma_rdec_interp(a_equality, m, C, T_ncdm0, reldec_factor)<H(a_equality)):
        return True 
    elif (Gamma_rdec_interp(a_equality, m, C, T_ncdm0, reldec_factor)>=H(a_equality)):
        return False

N_points_a = 3000
a_start = 1e-15
a_end = 1e0

def get_interpolator(mass,model,constant,T_ncdm0):

    #get factors
    reldec_factor, nrdec_factor = get_normalization_factors(T_ncdm0, mass)

    m = mass * 1e6 # in eV

    #assign blank functions depending on model
    if (model=='--C0'):
        Gamma_rdec_interp = Gamma_0_rdec_interp
        Gamma_nrdec_interp = Gamma_0_nrdec_interp
    elif (model == '--Cm'):
        Gamma_rdec_interp = Gamma_m_rdec_interp
        Gamma_nrdec_interp = Gamma_m_nrdec_interp
    elif (model == '--CV'):
        Gamma_rdec_interp = Gamma_V_rdec_interp
        Gamma_nrdec_interp = Gamma_V_nrdec_interp

    #calculate interaction rates (based on notebooks on NRdec! Basically carbon copies, except for updated T calculations...)
    a = np.logspace(np.log10(a_start),np.log10(a_end),N_points_a)
    Gamma_rel = Gamma_rdec_interp(a, m, constant, T_ncdm0, reldec_factor)
    Gamma_nrel = Gamma_nrdec_interp(a, m, constant, T_ncdm0, reldec_factor)

    #decide if reldec or nrdec is selected, return in proper units
    temperature_eV = T_ncdm_reldec(a, m, T_ncdm0)
    output_temp = unitConversion_eV_2_K(temperature_eV)
    if is_reldec(Gamma_rdec_interp, Gamma_nrdec_interp, m, constant, T_ncdm0, reldec_factor):
        output_reltime = unitConversion_time_evm1_2_Mpc(np.power(Gamma_rel, -1))
    else:
        output_reltime = unitConversion_time_evm1_2_Mpc(np.power(Gamma_nrel, -1))

    #print(output_temp)
    #print(output_reltime)

    return output_temp, output_reltime

def create_output_directories():
    path = os.getcwd() + '/reltime'
    if not os.path.exists(path):
        os.makedirs(path)

def output_tables(output_temp, output_reltime, mass, model, constant, psd_type, psd_file):

    #create table path
    output_path = os.getcwd() + '/reltime'
    output_file_prefix = '/relaxation_time_table_'
    file_number=0
    while os.path.exists(output_path+output_file_prefix+str(file_number)+'.dat'):
        file_number = file_number +1
    output_file_name = output_path+output_file_prefix+str(file_number)+'.dat'

    #wrie outputs
    with open(output_file_name,'w') as f:
        for i in range(N_points_a):
            f.write("%E \t %E \n" % (output_temp[i], output_reltime[i]))

    #create parameter file path. If it already exists, delete it
    output_file_name_params = output_path+output_file_prefix+str(file_number)+'_params.txt'
    if os.path.exists(output_file_name_params):
        os.remove(output_file_name_params)

    #write parameters
    m = mass * 1e6 #eV
    with open(output_file_name_params,'w') as f:
        f.write("mass= %E eV \n" % m)
        f.write("model= " + model[2:] + " \n")
        f.write("constant= %E \n" % constant)
        f.write("psd_type= " + psd_type[2:] + " \n")
        f.write("psd_file= " + psd_file + " \n")

    return output_file_name

def main():

    sys.stderr.write("Initialized interface_reltime.py. This is a console application of the routines in NRdec notebooks, accounting for T_ncdm0 != T_CMB. \n")

    #parse model, mass
    params = sys.argv[:]
    script = sys.argv[0]
    try:
        assert len(params)>=5
    except:
        return error_handler(-1)
    
    mass = float(params[1])
    constant = float(params[3])

    try:
        assert not(mass<0 or math.isnan(mass) or math.isinf(mass))
        assert not(constant<0 or math.isnan(constant) or math.isinf(constant))
    except:
        return error_handler(-2)

    model = params[2]
    try:
        assert model in ['--C0', '--CV', '--Cm']
    except:
        return error_handler(-3)
    if model in ['--CV', '--Cm']:
        constant = constant * (mass*1e6)**4

    psd_type = params[4]
    try:
        assert psd_type in ['--nonresonant', '--steriledm', '--thermal']
    except:
        return error_handler(-4)

    if (psd_type == '--steriledm'):
        try: 
            psd_file = params[5]
        except:
            return error_handler(-5)
    else:
        psd_file = ''

    sys.stderr.write('Parsed correctly all input parameters. Mass = ' + str(mass) + ' MeV, model = ' + model + ', coupling constant= ' + \
        str(constant) + ', psd model = ' + psd_type + ' psd_file = ' + psd_file + '. \n')

    #T_ncdm0 = get_Tncdm0_reldec(psd_type, os.getcwd() + '/' + psd_file, mass)
    T_ncdm0 = get_Tncdm0_reldec(psd_type, psd_file, mass)
    
    sys.stderr.write('Calculated T_ncdm0 = ' + str(T_ncdm0) + ' eV.\n')

    output_temp, output_reltime = get_interpolator(mass, model, constant, T_ncdm0)

    sys.stderr.write('Calculated relaxation time tables. \n')

    create_output_directories()

    output_file_name = output_tables(output_temp, output_reltime, mass, model, constant, psd_type, psd_file)

    sys.stderr.write('Tables written on file' + output_file_name + '. \n')

    return_function(0, output_file_name)

    return 0

if __name__ == '__main__':
   main()