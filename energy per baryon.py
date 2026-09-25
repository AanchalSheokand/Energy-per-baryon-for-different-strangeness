# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 17:09:02 2026

@author: aanch
"""

# -*- coding: utf-8 -*-
"""
Plot Energy per baryon vs chemical potential
for different strangeness fractions fs
with zoomed inset
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log, pi
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
# === Constants ===
g = 6.5
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
n = 0

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

    E_u = Eu
    E_d = Ed
    E_s = Es

    F1 = ((lambda_1+lambda_2)/4)*(su**3) \
        + ((2*m + lambda_1*(sd**2) + 2*lambda_1*(ss**2))/4)*su \
        - (c*sd*ss)/(2*sqrt(2)) - hud/2 - h3/2 \
        + (3/4)*((m_u*g)/(pi**2))*((p_u*E_u)-(m_u**2)*log((p_u+E_u)/m_u))

    F2 = ((lambda_1+lambda_2)/4)*(sd**3) \
        + ((2*m + lambda_1*(su**2) + 2*lambda_1*(ss**2))/4)*sd \
        - (c*su*ss)/(2*sqrt(2)) - hud/2 + h3/2 \
        + (3/4)*((m_d*g)/(pi**2))*((p_d*E_d)-(m_d**2)*log((p_d+E_d)/m_d))

    F3 = (lambda_1+lambda_2)*(ss**3) \
        + ((2*m + lambda_1*(su**2) + lambda_1*(sd**2))/2)*ss \
        - (c*su*sd)/(2*sqrt(2)) - hs \
        + (3/(2*sqrt(2)))*((m_s*g)/(pi**2))*((p_s*E_s)-(m_s**2)*log((p_s+E_s)/m_s))

    F4 = -(mw**2)*(omega)
    F5 = -(mw**2)*(ro)
    F6 = -(mp**2)*(phi)

    return [F1, F2, F3, F4, F5, F6]


# ================= PARAMETERS =================

s_values = [0,1,2,3]

fig, ax = plt.subplots(figsize=(15,20))

# Zoom inset
#axins = zoomed_inset_axes(ax, 3, loc='lower right')
axins = inset_axes(ax, width="65%", height="30%", loc="lower right",bbox_to_anchor=(0.5,0.5,2.5,2.5),
    bbox_transform=ax.transAxes)
# ================= MAIN LOOP =================

for s in s_values:

    mu_list=[]
    SA_val=[]
    rhobb=[]
    SA1_val=[]
    P_values=[]
    for i in np.linspace(0,25000,5000):
        A=i*0.0001
        rhob = i * 0.0001 * (fm**3)

        try:

            x=fsolve(fun,(90,90,90,0,0,0),args=(rhob,s,n))
            su,sd,ss,omega,ro,phi=x

            rho_u=(rhob/2)*(3-s-3*n-s*n)
            rho_d=(rhob/2)*(3-s+3*n+s*n)
            rho_s=s*rhob

            p_u=((pi**2)*rho_u)**(1/3)
            p_d=((pi**2)*rho_d)**(1/3)
            p_s=((pi**2)*rho_s)**(1/3)

            m_u=g*su/2
            m_d=g*sd/2
            m_s=g*ss/sqrt(2)

            Eu=sqrt(p_u**2+m_u**2)
            Ed=sqrt(p_d**2+m_d**2)
            Es=sqrt(p_s**2+m_s**2)

            mu=(Eu+Ed+Es)/3

            P1=(m/4)*(su**2+sd**2+2*(ss**2)) \
                -(c*su*sd*ss)/(2*sqrt(2)) \
                +(lambda_1/16)*((su**2+sd**2+2*(ss**2))**2) \
                +(lambda_2/16)*(su**4+sd**4+4*(ss**4)) \
                -hud*((su+sd)/2)-h3*((su-sd)/2)-hs*ss
            
            P=(1/(8*pi**2))*(
                2*p_u**3*Eu-3*m_u**2*p_u*Eu+3*m_u**4*log((p_u+Eu)/m_u)
            )+(1/(8*pi**2))*(
                2*p_d**3*Ed-3*m_d**2*p_d*Ed+3*m_d**4*log((p_d+Ed)/m_d)
            )+(1/(8*pi**2))*(
                2*p_s**3*Es-3*m_s**2*p_s*Es+3*m_s**4*log((p_s+Es)/m_s)
            )-P1-464*(fm**3)
            ED1= (3/(8*pi**2))*(
                2*p_u**3*Eu+ m_u**2*p_u*Eu-m_u**4*log((p_u+Eu)/m_u)
            )+(3/(8*pi**2))*(
                2*p_d**3*Ed+ m_d**2*p_d*Ed-m_d**4*log((p_d+Ed)/m_d)
            )+(3/(8*pi**2))*(
                2*p_s**3*Es+ m_s**2*p_s*Es-m_s**4*log((p_s+Es)/m_s)
            ) +P1+464*(fm**3)

            ED=-P+Eu*rho_u+Ed*rho_d+Es*rho_s

            SA= abs((ED)/(rhob))
            SA1=ED1/rhob
            mu_list.append(mu)
            SA_val.append(SA)
            SA1_val.append(SA1)
            rhobb.append((A))
            P_values.append(P/fm**3)
        except:
            continue

    ax.plot(mu_list,P_values,label=rf"$\eta={s}$") 
    #ax.plot(mu_list,P_values,label=rf"$\eta={s}$")
    axins.plot(mu_list,SA_val)
    ax.tick_params(axis='both',which='major', labelsize=30)
# ================= MAIN PLOT =================
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax.set_xlabel(r"$\mu_q$(MeV)",fontsize=30)
ax.set_ylabel(r"$E/A$ (MeV)",fontsize=30)
ax.set_xlim(330,460)# xmin, xmax
ax.set_ylim(0, 300)
ax.grid(True)
ax.legend(fontsize=30)

# ================= ZOOM REGION =================
axins.set_xlabel(r"$\mu_q(MeV)$", fontsize=30)
axins.set_ylabel(r"$E/A$ (MeV)",fontsize=30)
axins.set_xlim(330,360)
axins.set_ylim(800,1000)
ax.axhline(y=930.4, linestyle='--', label='Stability limit (930 MeV)')
axins.axhline(y=930.4, linestyle='--', label='Stability limit (930 MeV)')
axins.grid(True)

mark_inset(ax,axins,loc1=2,loc=4,fc="none",ec="0.5")

plt.title("Energy per baryon vs chemical potential")

plt.tight_layout()
plt.show()