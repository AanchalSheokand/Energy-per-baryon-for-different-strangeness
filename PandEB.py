# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:50:36 2026

@author: aanch
"""

# -*- coding: utf-8 -*-
"""
Energy per baryon and Pressure vs chemical potential
with zoomed inset for stability region
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log, pi
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# ================= CONSTANTS =================
g = 6
c = 4807.84
h_l = 1759728.23
h_s = 38072087.59
m = -1 * (306.26 ** 2)

lambda_1 = 13.49
lambda_2 = 46.48

fm = 197.3269804

hud = (120.73) ** 3
h3 = -1 * (78.31) ** 3
hs = (336.41) ** 3

mw = 782
mp = 1020

gv = 0
n = 0.1

# ================= GAP EQUATIONS =================
def fun(x, rhob, s, n):

    su, sd, ss, omega, ro, phi = x

    rho_u = (rhob / 2) * (3 - s - 3*n - s*n)
    rho_d = (rhob / 2) * (3 - s + 3*n + s*n)
    rho_s = s * rhob

    p_u = ((pi**2)*rho_u)**(1/3)
    p_d = ((pi**2)*rho_d)**(1/3)
    p_s = ((pi**2)*rho_s)**(1/3)

    m_u = g*su/2
    m_d = g*sd/2
    m_s = g*ss/sqrt(2)

    Eu = sqrt(p_u**2 + m_u**2)
    Ed = sqrt(p_d**2 + m_d**2)
    Es = sqrt(p_s**2 + m_s**2)

    F1 = ((lambda_1+lambda_2)/4)*(su**3) \
        + ((2*m + lambda_1*(sd**2) + 2*lambda_1*(ss**2))/4)*su \
        - (c*sd*ss)/(2*sqrt(2)) - hud/2 - h3/2 \
        + (3/4)*((m_u*g)/(pi**2))*((p_u*Eu)-(m_u**2)*log((p_u+Eu)/m_u))

    F2 = ((lambda_1+lambda_2)/4)*(sd**3) \
        + ((2*m + lambda_1*(su**2) + 2*lambda_1*(ss**2))/4)*sd \
        - (c*su*ss)/(2*sqrt(2)) - hud/2 + h3/2 \
        + (3/4)*((m_d*g)/(pi**2))*((p_d*Ed)-(m_d**2)*log((p_d+Ed)/m_d))

    F3 = (lambda_1+lambda_2)*(ss**3) \
        + ((2*m + lambda_1*(su**2) + lambda_1*(sd**2))/2)*ss \
        - (c*su*sd)/(2*sqrt(2)) - hs \
        + (3/(2*sqrt(2)))*((m_s*g)/(pi**2))*((p_s*Es)-(m_s**2)*log((p_s+Es)/m_s))

    F4 = -(mw**2)*(omega)
    F5 = -(mw**2)*(ro)
    F6 = -(mp**2)*(phi)

    return [F1, F2, F3, F4, F5, F6]

# ================= PARAMETERS =================
s_values = [0,0.1,0.2,1]

# ================= FIGURE =================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

# Zoom inset (for E/A)
axins = inset_axes(ax2, width="30%", height="28%", loc="upper right")

# ================= MAIN LOOP =================
for s in s_values:

    mu_list = []
    SA_val = []
    P_values = []
    su_val=[]
    sd_val=[]
    ss_val=[]
    P1_val=[]
    rhobb=[]
    for i in np.linspace(0, 30000,10000):
        A=i*(0.0001)
        rhob = i * 0.0001 * (fm**3)

        try:
            x = fsolve(fun, (90, 90, 90, 0, 0, 0), args=(rhob, s, n))
            su, sd, ss, omega, ro, phi = x

            rho_u = (rhob/2)*(3-s-3*n-s*n)
            rho_d = (rhob/2)*(3-s+3*n+s*n)
            rho_s = s*rhob

            p_u = ((pi**2)*rho_u)**(1/3)
            p_d = ((pi**2)*rho_d)**(1/3)
            p_s = ((pi**2)*rho_s)**(1/3)

            m_u = g*su/2
            m_d = g*sd/2
            m_s = g*ss/sqrt(2)

            Eu = sqrt(p_u**2+m_u**2)
            Ed = sqrt(p_d**2+m_d**2)
            Es = sqrt(p_s**2+m_s**2)

            mu = (Eu+Ed+Es)/3

            # ===== Pressure =====
            P1 = (m/4)*(su**2+sd**2+2*(ss**2)) \
                -(c*su*sd*ss)/(2*sqrt(2)) \
                +(lambda_1/16)*((su**2+sd**2+2*(ss**2))**2) \
                +(lambda_2/16)*(su**4+sd**4+4*(ss**4)) \
                -hud*((su+sd)/2)-h3*((su-sd)/2)-hs*ss

            P = (1/(8*pi**2))*(
                2*p_u**3*Eu-3*m_u**2*p_u*Eu+3*m_u**4*log((p_u+Eu)/m_u)
            ) + (1/(8*pi**2))*(
                2*p_d**3*Ed-3*m_d**2*p_d*Ed+3*m_d**4*log((p_d+Ed)/m_d)
            ) + (1/(8*pi**2))*(
                2*p_s**3*Es-3*m_s**2*p_s*Es+3*m_s**4*log((p_s+Es)/m_s)
            ) - P1 - 464.88250916189565*(fm**3)

            # ===== Energy Density =====
            ED = -P + Eu*rho_u + Ed*rho_d + Es*rho_s

            # ===== Energy per baryon =====
            SA = abs(ED/rhob)

            mu_list.append(mu)
            SA_val.append(SA)
            P_values.append(P/fm**3)
            su_val.append(su)
            sd_val.append(sd)
            ss_val.append(ss)
            rhobb.append(A/0.15)
            P1_val.append(-P1/fm**3)

        except:
            continue

    # ===== Plot =====
    ax1.plot(mu_list, P_values, label=rf"$\eta={s}$")
    ax2.plot(mu_list, SA_val, label=rf"$\eta={s}$")
    axins.plot(mu_list, SA_val)

# ================= PRESSURE PLOT =================
ax1.set_xlabel(r"$\mu_q$ (MeV)", fontsize=14)
ax1.set_ylabel(r"Pressure (MeV/fm$^3$)", fontsize=14)
ax1.set_title("Pressure vs Chemical Potential")
ax1.set_xlim(300,460)
ax1.set_ylim(0, 300)
ax1.grid(True)
ax1.legend()

# ================= ENERGY PER BARYON =================
ax2.set_xlabel(r"$\mu_q$ (MeV)", fontsize=14)
ax2.set_ylabel(r"$E/A$ (MeV)", fontsize=14)
ax2.set_title("Energy per Baryon vs Chemical Potential")
ax2.set_xlim(300,460)
ax2.set_ylim(700, 2000)
ax2.grid(True)
ax2.legend()

# Stability line (Bodmer-Witten condition)
ax2.axhline(y=930.4, linestyle='--')

# ================= ZOOM =================
axins.set_xlim(300, 360)
axins.set_ylim(800, 1000)
axins.axhline(y=930.4, linestyle='--')
axins.grid(True)

#mark_inset(ax2, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.tight_layout()
plt.show()