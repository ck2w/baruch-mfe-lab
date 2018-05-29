import numpy as np
from rough_bergomi_model import RoughBergomiModel
from local_volatility_model import LocalVolatilityModel
import math

def run_vol_of_vol(T, nuList, rho, H, fvc, r, S0, N):
    
    num_moneyness_rbergomi = 31
    num_moneyness_localvol = 31
    
    # K list of rBergomi
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_rbergomi)
    Klist_rBergomi = S0*np.exp(LogMoneyness_list)
    
    # K list of local vol model
    localVolN = N
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_localvol)
    Klist_localvol = S0*np.exp(LogMoneyness_list)
    
    euput_rBergomi_results = []
    amput_rBergomi_results = []
    euput_localvol_results = []
    amput_localvol_results = []
    
    for nu in nuList:
        
        # rBergomi
        rBergomi = RoughBergomiModel(H, rho, nu, fvc, r)
        Xmin, Xmax = rBergomi.create_price_tree(N, T, S0, EVOLVE_METHOD="dLogS", DEBUG=False)
        print("LogX span:", np.log(Xmin), np.log(Xmax))
        
        euput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=False, r=0.05)
        amput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=True, r=0.05)
        
        # Local vol
        localVol = LocalVolatilityModel(r=0.0)
        localVol.create_tree(localVolN, T, S0, rBergomi, N, False)

        euput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=False, r=0.05)
        amput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=True, r=0.05)
        
        euput_rBergomi_results.append(euput_rBergomi)
        amput_rBergomi_results.append(amput_rBergomi)
        euput_localvol_results.append(euput_localvol)
        amput_localvol_results.append(amput_localvol)
        
        print("*********************************************************************")
        
    return ( Klist_rBergomi, euput_rBergomi_results, amput_rBergomi_results, 
             Klist_localvol, euput_localvol_results, amput_localvol_results )

def run_eta(T, etaList, rho, H, fvc, r, S0, N):
    pass


def run_correlation(T, rhoList, nu, H, fvc, r, S0, N):
    
    num_moneyness_rbergomi = 31
    num_moneyness_localvol = 31
    
    # K list of rBergomi
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_rbergomi)
    Klist_rBergomi = S0*np.exp(LogMoneyness_list)
    
    # K list of local vol model
    localVolN = N
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_localvol)
    Klist_localvol = S0*np.exp(LogMoneyness_list)
    
    euput_rBergomi_results = []
    amput_rBergomi_results = []
    euput_localvol_results = []
    amput_localvol_results = []
    
    for rho in rhoList:
        
        # rBergomi
        rBergomi = RoughBergomiModel(H, rho, nu, fvc, r)
        Xmin, Xmax = rBergomi.create_price_tree(N, T, S0, EVOLVE_METHOD="dLogS", DEBUG=False)
        print("LogX span:", np.log(Xmin), np.log(Xmax))
        
        euput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=False, r=0.05)
        amput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=True, r=0.05)
        
        # Local vol
        localVol = LocalVolatilityModel(r=0.0)
        localVol.create_tree(localVolN, T, S0, rBergomi, N, False)

        euput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=False, r=0.05)
        amput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=True, r=0.05)
        
        euput_rBergomi_results.append(euput_rBergomi)
        amput_rBergomi_results.append(amput_rBergomi)
        euput_localvol_results.append(euput_localvol)
        amput_localvol_results.append(amput_localvol)
        
        print("*********************************************************************")
        
    return ( Klist_rBergomi, euput_rBergomi_results, amput_rBergomi_results, 
             Klist_localvol, euput_localvol_results, amput_localvol_results )
    
    
def run_roughness(T, Hlist, nu, rho, fvc, r, S0, N):
    
    num_moneyness_rbergomi = 31
    num_moneyness_localvol = 31
    
    # K list of rBergomi
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_rbergomi)
    Klist_rBergomi = S0*np.exp(LogMoneyness_list)
    
    # K list of local vol model
    localVolN = N
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_localvol)
    Klist_localvol = S0*np.exp(LogMoneyness_list)
    
    euput_rBergomi_results = []
    amput_rBergomi_results = []
    euput_localvol_results = []
    amput_localvol_results = []
    
    for H in Hlist:
        
        # rBergomi
        rBergomi = RoughBergomiModel(H, rho, nu, fvc, r)
        Xmin, Xmax = rBergomi.create_price_tree(N, T, S0, EVOLVE_METHOD="dLogS", DEBUG=False)
        print("LogX span:", np.log(Xmin), np.log(Xmax))
        
        euput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=False, r=0.05)
        amput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=True, r=0.05)
        
        # Local vol
        localVol = LocalVolatilityModel(r=0.0)
        localVol.create_tree(localVolN, T, S0, rBergomi, N, False)

        euput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=False, r=0.05)
        amput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=True, r=0.05)
        
        euput_rBergomi_results.append(euput_rBergomi)
        amput_rBergomi_results.append(amput_rBergomi)
        euput_localvol_results.append(euput_localvol)
        amput_localvol_results.append(amput_localvol)
        
        print("*********************************************************************")
        
    return ( Klist_rBergomi, euput_rBergomi_results, amput_rBergomi_results, 
             Klist_localvol, euput_localvol_results, amput_localvol_results )
    

def run_roughness_fixed_eta(T, Hlist, eta, rho, fvc, r, S0, N):
    
    num_moneyness_rbergomi = 31
    num_moneyness_localvol = 31
    
    # K list of rBergomi
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_rbergomi)
    Klist_rBergomi = S0*np.exp(LogMoneyness_list)
    
    # K list of local vol model
    localVolN = N
    LogMoneyness_list = np.linspace(-0.3, 0.3, num=num_moneyness_localvol)
    Klist_localvol = S0*np.exp(LogMoneyness_list)
    
    euput_rBergomi_results = []
    amput_rBergomi_results = []
    euput_localvol_results = []
    amput_localvol_results = []
    
    for H in Hlist:
        
        CH = math.sqrt((2*H*math.gamma(1.5-H))/(math.gamma(0.5+H)*math.gamma(2-2*H)))
        nu = eta/(2.0*CH)
        
        # rBergomi
        rBergomi = RoughBergomiModel(H, rho, nu, fvc, r)
        Xmin, Xmax = rBergomi.create_price_tree(N, T, S0, EVOLVE_METHOD="dLogS", DEBUG=False)
        print("LogX span:", np.log(Xmin), np.log(Xmax))
        
        euput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=False, r=0.05)
        amput_rBergomi,_ = rBergomi.price_option(N, T, Klist_rBergomi, opttype=-1, earlyexercise=True, r=0.05)
        
        # Local vol
        localVol = LocalVolatilityModel(r=0.0)
        localVol.create_tree(localVolN, T, S0, rBergomi, N, False)

        euput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=False, r=0.05)
        amput_localvol,_ = localVol.price_option(N, T, Klist_localvol, opttype=-1, earlyexercise=True, r=0.05)
        
        euput_rBergomi_results.append(euput_rBergomi)
        amput_rBergomi_results.append(amput_rBergomi)
        euput_localvol_results.append(euput_localvol)
        amput_localvol_results.append(amput_localvol)
        
        print("*********************************************************************")
        
    return ( Klist_rBergomi, euput_rBergomi_results, amput_rBergomi_results, 
             Klist_localvol, euput_localvol_results, amput_localvol_results )
    
    