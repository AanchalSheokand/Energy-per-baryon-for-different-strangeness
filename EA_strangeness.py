# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:44:20 2026

@author: aanch
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sqrt, log, pi
import math
# ================= PARAMETERS =================
def parameters():
    g=6.5
    pi=3.14
    c_1=4807.84
    h_l=1759728.23
    h_s=38072087.59
    m =-93795.1876
    lambda_1=13.49
    lambda_2=46.48
    fm=197.5
    return(g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)

(g,pi,c_1,h_l,h_s,m,lambda_1,lambda_2,fm)=parameters()
positive_infinity=math.inf

# ================= FUNCTIONS =================
def mesonic_l(sigma_l,sigma_s):
    return ((2*lambda_1+lambda_2)/2)*(sigma_l**3) + \
           (m-((sigma_s*c_1)/sqrt(2))+lambda_1*(sigma_s**2))*sigma_l - h_l

def mesonics__(sigma_l,sigma_s):
    return (lambda_1+lambda_2)*(sigma_s**3) + \
           (m+lambda_1*(sigma_l**2))*sigma_s - h_s - \
           ((c_1*(sigma_l**2))/(2*sqrt(2)))

def tempz_l(m_l,p_l,E_l):
    return ((3/2)*((m_l*g)/pi**2)*((p_l*E_l)-(m_l**2)*log((p_l+E_l)/m_l)))

def tempz_s(m_s,p_s,E_s):
    return ((3/(2*sqrt(2)))*((m_s*g)/pi**2)*((p_s*E_s)-(m_s**2)*log((p_s+E_s)/m_s)))

# ================= MAIN =================
fig, axes = plt.subplots(1,2, figsize=(20,10))
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#axins = inset_axes(axes[1], width="35%", height="35%", loc='lower left', bbox_to_anchor=(0, 0.1, 1, 1),  # shift upward (increase 0.1 → more up)
   # bbox_transform=axes[1].transAxes)
s_values =[ 3]

for s in s_values:

    mu_list = []
    P_list = []
    EB_list = []
    sl_list=[]
    ss_list=[]
    rhobb=[]
    P1_list=[]
    El_val=[]
    Es_val=[]
    ms_val=[]
    ml_val=[]
    def fun(x,rhob):
        sigma_l, sigma_s = x

        rho_s = s*rhob
        rho_l = ((3-s)/2)*rhob

        p_l=((pi**2)*rho_l)**(1/3)
        p_s=((pi**2)*rho_s)**(1/3)

        m_s=(g*sigma_s)/sqrt(2)
        m_l=(g*sigma_l)/2

        E_s=sqrt(p_s**2+m_s**2)
        E_l=sqrt(p_l**2+m_l**2)

        F1=mesonic_l(sigma_l,sigma_s)+tempz_l(m_l,p_l,E_l)
        F2=mesonics__(sigma_l,sigma_s)+tempz_s(m_s,p_s,E_s)

        return [F1,F2]

    for i in range(1,300):
        A=i*0.01
        rhob=A*(fm**3)

        rho_s = s*rhob
        rho_l = ((3-s)/2)*rhob

        p_l=((pi**2)*rho_l)**(1/3)
        p_s=((pi**2)*rho_s)**(1/3)

        sl,ss=fsolve(fun,(90,90),(rhob))

        m_l=g*sl/2
        m_s=g*ss/sqrt(2)

        E_s=sqrt(p_s**2+m_s**2)
        E_l=sqrt(p_l**2+m_l**2)

        # ===== Chemical potential =====
        mu_B = (2*E_l + E_s)/3

        # ===== Pressure =====
        P1 = (m/2)*(sl**2+(ss**2)) \
            -(c_1*sl*sl*ss)/(2*sqrt(2)) \
            +(lambda_1*(sl**2)*(ss**2))/(2)\
            +((2*lambda_1+lambda_2)*sl**4/(8))\
            +((lambda_1+lambda_2)*(ss**4))/(4)-h_l*sl-h_s*ss

        P = (2/(8*pi**2))*(
            2*p_l**3*E_l-3*m_l**2*p_l*E_l+3*m_l**4*log((p_l+E_l)/m_l)
        ) + (1/(8*pi**2))*(
            2*p_s**3*E_s-3*m_s**2*p_s*E_s+3*m_s**4*log((p_s+E_s)/m_s)
        ) - P1 - 463.64765556922464*(fm**3)

        # ===== Energy density =====
        E1 = -P + 2*rho_l*E_l + rho_s*E_s

        EB = abs(E1 / rhob)
        sl_list.append(sl)
        ss_list.append(ss)
        mu_list.append(mu_B)
        P_list.append(P/fm**3)
        EB_list.append(EB)
        rhobb.append(A)
        P1_list.append(E1/fm**3)
        El_val.append(E_l)
        Es_val.append(E_s)
        ms_val.append(m_s)
        ml_val.append(m_l)
        eos_array = np.column_stack((rhobb,P_list, P1_list))
        np.savetxt("EAS.dat", eos_array)
    # ===== Plot =====
    axes[0].plot(mu_list, EB_list, label=f"$\eta$={s}")
    axes[1].plot(ml_val, EB_list, label=f"$\eta$={s}")
    #axins.plot(ml_val, EB_list)
    #axins.tick_params(axis='both',which='major', labelsize=30) 
    axes[0].tick_params(axis='both',which='major', labelsize=30)
    axes[1].tick_params(axis='both',which='major', labelsize=30)
# ================= LABELS =================
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

plt.yticks(fontsize=30)
axes[0].set_xlabel(r"$m_s(MeV)$",fontsize=30)
axes[0].set_ylabel(r"E/A(MeV)", fontsize=30)
#axes[0].set_title("Pressure vs Chemical Potential")
axes[0].grid(True)
axes[0].set_xlim(330,460)
axes[0].set_ylim(900, 1400)
axes[1].set_xlabel(r"$m_l(MeV)$ ", fontsize=30)
axes[1].set_ylabel(r"$E/A (MeV)$", fontsize=30)
#axes[1].set_title("Energy per Baryon vs Chemical Potential")
axes[1].grid(True)
axes[1].set_xlim(0,310)
axes[1].set_ylim(900, 1400)

# ================= ZOOM INSET =================

#axins.set_xlim(0,310)
#axins.set_ylim(800,1000)
#axins.grid(True)
axes[1].axhline(y=930.4, linestyle='--')
#axins.axhline(y=930.4, linestyle='--')
axes[0].axhline(y=930.4, linestyle='--')

# ================= LEGEND =================
axes[0].legend(fontsize=30)
#axes[1].legend(fontsize=20)

plt.tight_layout()
plt.show()