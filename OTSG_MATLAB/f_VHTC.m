function [UA,m_dot_Pri,m_dot_Sec] = f_VHTC(PD_Ratio,T_Pri_m,P_Sec,T_Sec_m,T_Sec_feed)
%Give a dummi Valve
%x_f = 0;
%For the Primary side:
D_Tube = 0.0196;
Delta_Tube = 0.0012;
m_dot_Pri = 1000;
CSA_SG = 2.751;
P_Pri = 145;
Power = 2e8;
A_Pri = (PD_Ratio^2-0.25*pi)/PD_Ratio^2*CSA_SG;
D_h_Pri = 0.0253;
h_Pri = XSteam('h_pT',P_Pri,T_Pri_m);
pho_Pri = XSteam('rho_ph',P_Pri,h_Pri);
my_Pri = XSteam('my_ph',P_Pri,h_Pri);
Cp_Pri = XSteam('Cp_ph',P_Pri,h_Pri);
tc_Pri = XSteam('tc_ph',P_Pri,h_Pri);
V_Pri = m_dot_Pri / (pho_Pri*A_Pri);
Re_Pri = V_Pri * D_h_Pri * pho_Pri / my_Pri;
Pr_Pri = my_Pri * Cp_Pri * 1000 / tc_Pri;
Nu_Pri = f_Nu_Gnielinski(Re_Pri,Pr_Pri);
h_Pri = Nu_Pri * tc_Pri / D_h_Pri;
%For the Secondary side:
D_h_Sec = D_Tube;
A_Sec = (0.25*pi)/PD_Ratio^2*CSA_SG;
Num_Tube = fix(A_Sec / (0.25 * pi * D_Tube^2));
m_dot_Sec = Power / 1000 / (XSteam('h_pT',P_Sec,XSteam('Tsat_p',P_Sec)+30)-XSteam('h_pT',P_Sec,T_Sec_feed));
h_Sec = XSteam('h_pT',P_Sec,T_Sec_m);
pho_Sec = XSteam('rho_ph',P_Sec,h_Sec);
my_Sec = XSteam('my_ph',P_Sec,h_Sec);
Cp_Sec = XSteam('Cp_ph',P_Sec,h_Sec);
tc_Sec = XSteam('tc_ph',P_Sec,h_Sec);
V_Sec = m_dot_Sec / (pho_Sec*A_Sec);
Re_Sec = V_Sec * D_h_Sec * pho_Sec / my_Sec;
Pr_Sec = my_Sec * Cp_Sec * 1000 / tc_Sec;
Nu_Sec = f_Nu_Gnielinski(Re_Sec,Pr_Sec);
h_Sec = Nu_Sec * tc_Sec / D_h_Sec;
%For the Tube Wall:
tc_TubeWall = 18.3;
L = 0.1;
R_w = log((D_Tube+2*Delta_Tube)/D_Tube)/(2*pi*L*tc_TubeWall);
%The Overall Heat Transfer Coefficient:
UA = 1 / (R_w + 1/(h_Pri*(pi*(D_Tube+2*Delta_Tube)*L)) + 1/(h_Sec*(pi*D_Tube*L)));
