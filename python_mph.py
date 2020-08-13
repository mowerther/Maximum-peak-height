# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:26:45 2020

@author: Mortimer
"""

#####################################################################################################################
# Maximum-peak height (MPH) Sentinel-3 OLCI Python example implementation - Mortimer Werther (2020)                  #
# The code is based on the publication by M. Matthews and D. Odermatt 2015:                                         #
# "Improved algorithm for routine monitoring of cyanobacteria and eutrophication in inland and near-coastal waters" #
# https://doi.org/10.1016/j.rse.2014.10.010                                                                         #
#####################################################################################################################

# dependencies
import math
import numpy as np

# Splitting the calcuation of MPH CHL into two functions:
# 1. Required MPH run procedure parameters such as the peaks and positions
# 2. MPH run procedure using the parameters from 1.

def calc_mph_paras(rrs_620, rrs_665, rrs_681, rrs_709, rrs_753, rrs_885):
    
    '''
    Calculation of MPH parameters
    Requires satellite-derived BRR or in situ reflectance at 620nm, 665nm, 681nm, 709nm, 753nm and 885nm (OLCI)
    '''
    
    rmax_0 = max(rrs_681, rrs_709)
    
    if rrs_681 >= rrs_709:
        lambda_rmax_0 = 681
    else:
        lambda_rmax_0 = 709
    
    rmax_1 = max(rrs_681, rrs_709, rrs_753)
    
    if (rrs_681 >= rrs_709) & (rrs_681 >= rrs_753):
        lambda_rmax_1 = 681
    elif (rrs_709 >= rrs_681) & (rrs_709 >= rrs_753):
        lambda_rmax_1 = 709
    else:
        lambda_rmax_1 = 753
        
    # NDVI
    ndvi = (rrs_885 - rrs_665) / (rrs_885 + rrs_665)
    
    # SIPF, SICF, BAIR
    sipf = rrs_665 - rrs_620 - ((rrs_681 - rrs_620) * (665 - 620) / (681 - 620))
    sicf = rrs_681 - rrs_665 - ((rrs_709 - rrs_665) * (681 - 665) / (709 - 665))
    bair = rrs_709 - rrs_665 - ((rrs_885 - rrs_665) * (709 - 665) / (885 - 665))
    
    # MPH_0, MPH_1
    s
    mph_0 = rmax_0 - rrs_665 - ((rrs_885 - rrs_665) * (lambda_rmax_0 - 665) / (885 - 665))
    mph_1 = rmax_1 - rrs_665 - ((rrs_885 - rrs_665) * (lambda_rmax_1 - 665) / (885 - 665))

    return mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0

def calc_mph_chl(mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0, mph_floatthres, mph_cyanomax):
    
    """
    Calculation of MPH CHL
    Requires previously calculated parameters AND mph_floatthres and mph_cyanomax (high values, see paper)
    Note: obviously out-comment or delete print statements if you intent to use this code.
    """
    
    if lambda_rmax_1 == 753:
        print('Right side of if-condition.')
        # MPH >= 0.02 or NDVI >0.2
        if (mph_1 >= 0.02 or ndvi >= 0.2):
            float_flag = 1
            adj_flag = 0
            # SICF < 0 and SIPF > 0
            if (sicf < 0 and sipf > 0):
                cyano_flag=1
                print('Flag: floating cyanobacteria true')
                chl_mph = 22.44 * math.exp(35.79 * mph_1)
                print('CHL MPH is: ' + str(chl_mph))
                if chl_mph > mph_floatthres:
                    float_flag=1
                    print('Floating cyanobacteria')
                else:
                    print('Immersed cyanobacteria')
            # SICF >=0 or SIPF <=0 
            elif (sicf >= 0 or sipf <= 0):
                cyano_flag = 0
                chl_mph = np.nan
                print('Floating vegetation')
    
        # Continuation right side
        elif (mph_1 < 0.02 and ndvi < 0.2):
            float_flag = 0
            adj_flag = 1
            print('Flag: adjacent true')
            cyano_flag = 0
            print('Immersed eukaryotes')
    
            chl_mph = 5.24 * 10 ** 9 * mph_0 ** 4 - 1.95 * 10 ** 8 * mph_0 ** 3 + 2.46 * 10 ** 6 * mph_0 ** 2 + 4.02 * 10 ** 3 * mph_0 + 1.97
    
    # Left side of if-condition
    else:
        print('Left side of if-condition.')
        float_flag = 0
        adj_flag = 0
        
        # Left side of 2nd if-condition
        if (sicf >= 0 or sipf <= 0 or bair <= 0.002):
                print('Left 2nd if-condition')
                cyano_flag=0
                print('Immersed eukaryotes')
                chl_mph = 5.24 * 10 ** 9 * mph_0 ** 4 - 1.95 * 10 ** 8 * mph_0 ** 3 + 2.46 * 10 ** 6 * mph_0 ** 2 + 4.02 * 10 ** 3 * mph_0 + 1.97
        
        # Right side of 2nd if-condition
        elif (sicf <= 0 and sipf > 0 and bair > 0.002):
            print('Right 2nd if-condition')
            cyano_flag = 1
            print('Flag: cyanobacteria true')
            chl_mph = 22.44 * math.exp(35.79 * mph_1)
            if chl_mph > mph_floatthres:
                    float_flag = 1
                    print('Floating cyanobacteria')
                    if chl_mph > mph_cyanomax:
                        chl_mph = mph_cyanomax
                        print('MPH chl maximum reached.')

    return chl_mph, cyano_flag, float_flag, adj_flag


# Sample Rrs (Sentinel - 3 OLCI bands)
rrs_620 = 0.0014
rrs_665 = 0.0020
rrs_681 = 0.0042
rrs_709 = 0.0025
rrs_753 = 0.0010
rrs_885 = 0.0003

# Call parameter function
mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0 = calc_mph_paras(rrs_620, rrs_665, rrs_681, rrs_709, rrs_753, rrs_885)

# Call MPH run procedure
# Feel free to choose different floatthres and cyanomax values
chl_mph, adj_flag, cyano_flag, float_flag = calc_mph_chl(mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0, 200, 350)

print("MPH CHL is: " + str(chl_mph))
# should be ~22.3
print("Flags - Adj flag: " + str(adj_flag)+ ", "+ "Cyano Flag: " + str(cyano_flag) +", "+ "Float flag: " + str(float_flag))
# should be: Flags - Adj flag: 0, Cyano Flag: 0, Float flag: 0
