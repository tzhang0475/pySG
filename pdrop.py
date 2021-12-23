# module for solve pressure drop in one verticali uniform tube of single phase or two phase HEM model

# note: flow velocity vector vf: 1 for upwards flow and -1 for downwards flow

# lib to use
import numpy as np
import math
from iapws import IAPWS97
import fluid_data as fd
import hcoeff as hc

# define constant
g = 9.81


# single phase pressure drop calculation
def pd_single(mdot_n_ip1,mdot_n_i,Tdata,Pdata,rho_n_ip1,rho_n_i,x,dh,A_xsc,l_cell,vf):

#	calulate viscosity and density
	mu = IAPWS97(P = Pdata, T = Tdata).mu
	rho = IAPWS97(P = Pdata, T = Tdata).rho
#	calculate average mass flow
	mdot = 1/2 * (mdot_n_ip1 + mdot_n_i) 
#	calculate Re 
	Re = hc.Reynold(dh,mdot,A_xsc,mu)
#	calculate friction factor 
	f = hc.f_colebrook(Re)
#	calculate friction pressure drop
	dp_f = f* l_cell/dh * mdot**2/(2*rho*A_xsc**2) 
#	calculate the accelorate pressure drop
	dp_a = mdot_n_ip1**2/(A_xsc**2*rho_n_ip1) - mdot_n_i**2/(A_xsc**2*rho_n_i)
#	calculate the gravity pressure drop
	dp_g = 1/2*(rho_n_ip1+rho_n_i) * g * l_cell

#	totoal pressure drop
	result = (dp_f + dp_a + vf*dp_g)/1e6 #sum up pressure drop and convert to MPa

	return result

# single phase dynamic pressure drop
def pd_dyn_single(u_n_ip1,u_n_i,rho_n_ip1,rho_n_i):
#	result = mdot**2/(A_xsc**2*rho_n_ip1) - mdot**2/(A_xsc**2*rho_n_i)
#	result = mdot_n_ip1**2/(A_xsc**2*rho_n_ip1) - mdot_n_i**2/(A_xsc**2*rho_n_i)
	result = rho_n_ip1*u_n_ip1**2 - rho_n_i*u_n_i**2

	return result #in Pa

# calculate two phase flow  viscosity
def mu_tw(Pdata,xdata):
	# calculate viscosity of liquid and vapor
	mu_l = IAPWS97(P = Pdata, x = 0.0).mu
	mu_g = IAPWS97(P = Pdata, x = 1.0).mu

	# calculate ratio of  two phase viscosity and liquid viscosity , based on McAdans et la
	mu_tw = (mu_l*(1+xdata*(mu_l/mu_g - 1)))**(-1)
	# calculate ratio of  two phase viscosity and liquid viscosity , based on Cichitti et la
#	mu_tw = (1+xdata*(mu_g/mu_l - 1))*mu_l
	# calculate  two phase flow  viscosity, based on lin et la
#	mu_tw = mu_g*mu_l/(mu_g + xdata**1.4*(mu_l - mu_g))

# Pressure drop based on liquid alone, l
def dp_l(Pdata,mdot,l_cell,dh,A_xsc): 
#	calulate viscosity and density
	mu = IAPWS97(P = Pdata, x = 0).mu
	rho = IAPWS97(P = Pdata, x = 0).rho
#	calculate Re 
	Re = hc.Reynold(dh,mdot,A_xsc,mu)
#	calculate friction factor 
	f = hc.f_colebrook(Re)
#	calculate friction pressure drop
	dp_f = f* l_cell/dh * mdot**2/(2*rho*A_xsc**2) 

	result = dp_f

	return result
# Pressure drop based on vapor alone, g 
def dp_g(Pdata,mdot,l_cell,dh,A_xsc): 
#	calulate viscosity and density
	mu = IAPWS97(P = Pdata, x = 1).mu
	rho = IAPWS97(P = Pdata, x = 1).rho
#	calculate Re 
	Re = hc.Reynold(dh,mdot,A_xsc,mu)
#	calculate friction factor 
	f = hc.f_colebrook(Re)
#	calculate friction pressure drop
	dp_f = f* l_cell/dh * mdot**2/(2*rho*A_xsc**2)

	result = dp_f

	return result 
# Pressure drop based on liquid only, lo 
def dp_lo(Pdata,xdata,mdot,l_cell,dh,A_xsc): 
#	calulate viscosity and density
	mu = mu_tw(Pdata,xdata)
	rho = IAPWS97(P = Pdata, x = xdata).rho
#	calculate Re 
	Re = hc.Reynold(dh,mdot,A_xsc,mu)
#	calculate friction factor 
	f = hc.f_colebrook(Re)
#	calculate friction pressure drop
	dp_f = f* l_cell/dh * mdot**2/(2*rho*A_xsc**2) 

	result = dp_f

	return result
# Pressure drop based on vapor only, go
# homogenous model of two phase friction pressure drop
# Armand correlation for two phase friction pressure drop multiplier
def phi_arm(Pdata, xdata):
	# calculate density of vapor and two phase flow
	rho_g = IAPWS97(P = Pdata, x = 1.0).rho
	rho_tw = IAPWS97(P = Pdata, x = xdata).rho

	# calculate two phase flow void fraction
	alf = xdata * rho_tw/rho_g
	# two phase pressure drop multiplier by Armand corrolation, from cobra en
	if alf >= 0 and alf < 0.6:
		phi_arm = (1-xdata)**2/(1-alf)**1.42
	elif alf >= 0.6 and alf < 0.9:
		phi_arm = 0.478 * (1-xdata)**2/(1-alf)**2.2
	elif alf >= 0.9 and alf <= 1.0:
		phi_arm = 1.73*(1-xdata)**2/(1-alf)**1.64
	
	result = phi_arm

	return result
# Lockhart-Martinelli correlation for two phase friction pressure drop multiplieri, liquid alone 
def phi_l_LM(mdot, Pdata, xdata,dh,l_cell,A_xsc):
	# calculate Martinelli parameter
	phi_l_sqr = dp_l(Pdata,mdot,l_cell,dh,A_xsc)/l_cell
	phi_g_sqr = dp_g(Pdata,mdot,l_cell,dh,A_xsc)/l_cell
	chi_sqr = phi_g_sqr/phi_l_sqr

        # calulate viscosity of liquid and vapor
	mu_l = IAPWS97(P = Pdata, x = 0).mu
	mu_g = IAPWS97(P = Pdata, x = 1).mu
	# calculate mass flow for each phase
	mdot_l = (1-xdata) * mdot
	mdot_g = xdata * mdot
        # calculate Re 
	Re_l = hc.Reynold(dh,mdot_l,A_xsc,mu_l)
	Re_g = hc.Reynold(dh,mdot_g,A_xsc,mu_g)

	# calculate value of C
	if Re_l >= 2400 and Re_g >= 2400:
		C = 20
	elif Re_l < 2400 and Re_g >= 2400:
		C = 12
	elif Re_l >= 2400 and Re_g < 2400:
		C = 10
	elif Re_l < 2400 and Re_g < 2400:
		C = 5

	result = 1+ C/math.sqrt(chi_sqr) + 1/chi_sqr

	return result
# Chisholm correlation for two phase friction pressure drop multiplier
#def phi_lo_Chris(Pdata, xdata):
# HEM model two phase flow pressure drop
def pd_twHEM(mdot_n_ip1,mdot_n_i,Tdata,Pdata,rho_n_ip1,rho_n_i,x_i1,x_i0,dh,A_xsc,l_cell,vf):
#	mean vapor fraction and mass flow
	xdata = 1/2 * (x_i0 + x_i1)
	mdot = 1/2 * (mdot_n_ip1 + mdot_n_i)	

#       calculate density of liquid and vapor
	rho_l = IAPWS97(P = Pdata, x = 0.0).rho
	rho_g = IAPWS97(P = Pdata, x = 1.0).rho
	rho_tw = IAPWS97(P = Pdata, x = xdata).rho
#       calculate viscosity of liquid and vapor
	mu_l = IAPWS97(P = Pdata, x = 0.0).mu
	mu_g = IAPWS97(P = Pdata, x = 1.0).mu

#       calculate liquid and vapor part Re
	Re_l = hc.Reynold(dh,mdot,A_xsc,mu_l)
#       calculate liquid friction factor
	f_l = hc.f_colebrook(Re_l)

#	r_mu = mu_tw/mu_l
#	Re_tw = hc.Reynold(dh,mdot,A_xsc,mu_tw)

#	power parameter
#	if Re_tw < 2000:
#		n = 1
#		C = 0.316
#	elif Re_tw >= 2000 and Re_tw <= 20000:
#		n = -0.25
#		C = 0.316
#	elif Re_tw >= 20000:
#		n = -0.2
#		C = 0.184
#	calculate the two phase flow friction factor
#	f_tp = f_l * r_mu**(0.2)
#	f_tp = C * Re_tw**(n)
#	f_r = f_tp/f_l

	# two phase pressure drop multiplier by Armand corrolation, from cobra en
#	phi_arm = phi_arm(Pdata,xdata)
	# Lockhart-Martinelli correlation for two phase friction pressure drop multiplieri, liquid alone 
	phi_l = phi_l_LM(mdot, Pdata, xdata,dh,l_cell,A_xsc)
	# single liquid phase friction presure drop
	dp_fl = dp_l(Pdata,mdot,l_cell,dh,A_xsc)
#	calculate the friction pressure drop
	dp_f = phi_l * dp_fl
#	dp_f = phi_arm * f_l/dh * (mdot**2/(2*A_xsc*rho_l)) * l_cell
#	dp_f = f_tp/dh * (mdot**2/(2*A_xsc*rho_tw)) * l_cell * (1 + xdata * (rho_l/rho_tw))

#       calculate acceleration pressure drop
#	dp_a = (mdot/A_xsc)**2 *(1/rho_l + (1/rho_g - 1/rho_l)*xdata) * l_cell
	dp_a = mdot_n_ip1**2/(A_xsc**2*rho_n_ip1) - mdot_n_i**2/(A_xsc**2*rho_n_i)

#       calculate gravity preesure drop
	dp_g = IAPWS97(P = Pdata, x = xdata).rho * g * l_cell

#	total pressure drop
	result = (dp_f + dp_a + vf*dp_g)/1e6 #sum up pressure drop and convert to MPa

	return result

# two phase dynamic preesure drop
def pd_dyn_twHEM(mdot_n_ip1,mdot_n_i,Pdata,rho_ip1,rho_i,x_i1,x_i0,A_xsc):
	mdot = 1/2* (mdot_n_ip1 + mdot_n_i)
#       calculate acceleration pressure drop
#	result = (mdot/A_xsc)**2 *(1/rho_l + (1/rho_g - 1/rho_l)*xdata) * l_cell
	result = mdot_n_ip1**2/(A_xsc**2*rho_ip1) - mdot_n_i**2/(A_xsc**2*rho_i)
#	result = (mdot/A_xsc) * (1/rho_ip1 - 1/rho_i)

	return result #in Pa
