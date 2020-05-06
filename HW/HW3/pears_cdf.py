# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 08:56:11 2018

@author: Annalise
"""
import numpy as np
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import gamma
from scipy.stats import f

def pearson_fit(end, mu, sig, beta1, beta2):
    k=1
    plotf=1
    output=1
    method='G.Q.'

    if end[1] == np.inf:
        case = 1
    elif end[0] == -np.inf:
        case = 2
    else:
        case = 3
    end = [(i - mu) / sig for i in end] #normalize end
    a = (4. * beta2 - 3. * beta1) #/ (10. * beta2 - 12. * beta1 - 18.) * pow(sig,2)
    b = np.sqrt(mu * beta1) * (beta2 + 3.) #/ (10. * beta2 - 12. * beta1 - 18.)
    c = (2. * beta2 - 3. * beta1 - 6.) #/ (10. * beta2 - 12 * beta1 - 18.)
    print(a,b,c)
    #figure out pearson type
    if abs(b) < 1e-5: 
        if abs(beta2 - 3) < 1e-5:#and abs(c) < 1e-5 and a > 0.: #pearson 0 - Gaussian
            pears_type = 0
            a1 = 0.
            a2 = 0.
        elif beta2 < 3.:
            pears_type = 2
            a1 = -np.sqrt(abs(a/c))
            a2 = a1 * -1
        elif beta2 > 3:
            pears_type = 7
            a1 = -np.sqrt(abs(a/c))
            a2 = a1 * -1
    elif abs(c) < 1e-5:
        pears_type = 3
        a1 = a * -1
        a2 = a1 * 1
    else:
        kappa = pow(b,2) / (4.*a*c)
        if kappa < 0.:
            pears_type = 1
        elif kappa < 1 - 1e-10:
            pears_type = 4
        elif kappa <= 1 + 1e-10:
            pears_type = 5
        else:
            pears_type = 6
        a1 = (-b + np.sqrt(pow(b,2) - 4.*a*c)) / (2. * c)
        a2 = (-b - np.sqrt(pow(b,2) - 4.*a*c)) / (2. * c)
        #reorganize roots so that a1 < 0 < a2
        if a1 > 0. and a2 < 0.: #sort the roots so a1 < 0 < a2
            dummy = a1 * 1.
            a2 = a1 * 1.
            a1 = dummy * 1.
        elif a1 < 0. and a2 > 0.:
            pass
        else:
            return ('error in roots')
    #use type to determine cdf
    denom = (10. * beta2 - 12. * beta1 - 18.)
    if abs(denom) > np.sqrt(2.2251e-308): #don't devide by 0
        a = a / denom
        b = b / denom
        c = c / denom
        coefs = [a, b, c]
    else:
        pears_type = 1
        coefs = [np.inf, np.inf, np.inf]
    
    if method == 'MCS':
        pears_type = 8
    #if pearson normal distribution
    if pears_type == 0:
        m1 = 0.
        m2 = 1.
        p = norm.cdf(end[1],m1,m2) - norm.cdf(end[0],m1,m2)
        inv1 = norm.ppf(p)
        inv2 = norm.ppf(norm.cdf(end[1],m1,m2))

    elif pears_type == 1:
        if abs(denom) > np.sqrt(2.2251e-308): #don't devide by 0
            m1 = (b + a1) / (c * (a2 - a1))
            m2 = -(b + a2) / (c * (a2 - a1))
        else:
            m1 = b / (c * (a2 - a1))
            m2 = -b / (c * (a2 - a1))
        end = [(i - a1) / (a2 - a1) for i in end]
        p = beta.cdf(end[1],m1 + 1, m2 + 1) - beta.cdf(end[0],m1 + 1, m2 + 1)
    elif pears_type == 2:
        m1 = (b + a1) / (c * 2 * abs(a1))
        m2 = ma * 1.
        end = [(i - a1) / (2*abs(a1)) for i in end]
        p = beta.cdf(end[1],m1 + 1, m2 + 1) - beta.cdf(end[0],m1 + 1, m2 + 1)
    elif pears_type == 3:
        m1 = (a / b - b) / b
        m2 = m1 * 1.
        end = [(i - a1) / b for i in end]
        p = gamma.cdf(end[1],m1 + 1) - gamma.cdf(end[0],m1 + 1)
    elif pears_type == 4:
        end = [(i * sig + mu) for i in end]
        r = 6. * (beta2 - beta1 - 1.) / (2*beta2 - 3*beta1 - 6.)
        m1 = 1. + r / 2.
        m2 = -r*(r - 2.) * np.sqrt(beta1) / np.sqrt(16.*(r - 1.) - beta1*pow(r - 2.,2))
        a = np.sqrt(pow(sig,2) * (16*(r - 1.) - beta1*pow(r - 2.,2)))/4.
        lam = mu - ((r - 2.)*np.sqrt(beta1) * sig)/4.
        '''Python doesn't seem to have a built in function for pearson 4
        Our homework doesn't use perason 4
        so I'm just going to leave this out for now'''
        #if case == 1:
        #    p = 1 - 
        #elif case == 2:
        #    0 = 
        #else:
        #    p = 
    elif pears_type == 5:
        C1 = b / (2.*c)
        end = [-((b - C1)/c)/(i + C1) for i in end]
        m1 = c * 1.
        m2 = 0.
        p = gamma.cdf(end[1],1./c - 1) - gamma.cdf(end[0],1./c - 1)
    elif pears_type == 6:
        m1 = (a1 + b) / (c*(a2 - a1))
        m2 = -(a2 + b) / (c*(a2 - a1))
        if a2 < 0.:
            nu1 = 2.*(m2 + 1.)
            nu2 = -2*(m1 + m2 + 1.)
            end = [(i - a2)/(a2 - a1)*(nu2 / nu1) for i in end]
            p = f.cdf(end[1],nu1,nu2) - f.cdf(end[0],nu1,nu2)
        else:
            nu1 = 2.*(m1 + 1.)
            nu2 = -2.*(m1 + m2 + 1.)
            end = [(i - a1)/(a1 - a2)*(nu2 / nu1) for i in end]
            p = f.cdf(end[1],nu1,nu2) - f.cdf(end[0],nu1,nu2)
    elif pears_type == 7:
        m1 = 1./c - 1.
        end = [(i/np.sqrt(a/(1. - c))) for i in end]
        m2 = 0.
        p = t.cdf(end[1],m1) - t.cdf(end[0],m1)
    #print(m1, m2)
    print('distribution type: ',pears_type)
    return p, pears_type