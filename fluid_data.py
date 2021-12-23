# module used for generator fluid data calculate for heat transfer coefficient

# lib in use
import numpy as np
from iapws import IAPWS97 

def fdata(Tdata,Pdata,xdata): # T for fluid temperature in K, P for fluid pressure in MPa, x for steam fraction(0<=x<=1)
#	initial array and parameters
	result = np.ones(6)    
	rho = 0   
	k = 0
	mu = 0
	cp = 0
	h = 0
	x_tt= 0 

	   

	if (xdata>0)&(xdata<1.0):
#		gas density
		rho_g = IAPWS97(P = Pdata, x = 1).rho
#		gas viscosity
		mu_g  = IAPWS97(P = Pdata, x = 1).mu 
#		fluid density
		rho_f = IAPWS97(P = Pdata, x = 0).rho
#		fluid viscosity
		mu_f = IAPWS97(P = Pdata, x = 0).mu
#		two phase flow density
		rho = IAPWS97(P = Pdata, x= xdata).rho
	       
#		fluid conductivity in W/(K.m)
		k = 0 # no data for conductivity for two phase flow
	   
#		fluid viscosity in Pa.s
		mu = 0
	   
#		fluid heat capacity in J/(kg.K)
		cp = 0  # no heat capacity for two phase flow

#		steam enthalpy
		h = IAPWS97(P = Pdata, x = xdata).h
#		two phase correction factor	
	#	x_tt = ((1-xdata)/xdata)**(0.9) * (rho_g/rho_f)**0.5 * (mu_f/mu_g)**(0.1) 
		x_tt = 1/((xdata/(1-xdata))**(0.9) * (rho_f/rho_g)**(0.5) * (mu_g/mu_f)**(0.1))
	else:
#		liquid perproties
	
#		fluid density, in kg/m^3 
		rho = IAPWS97(P = Pdata, T = Tdata).rho
	       
#		fluid conductivity in W/(K.m)
		k = IAPWS97(P = Pdata, T = Tdata).k
	   
#		fluid viscosity in Pa.s
		mu = IAPWS97(P = Pdata, T = Tdata).mu
	   
#		fluid heat capacity in J/(kg.K)
		cp = IAPWS97(P = Pdata, T = Tdata).cp * 1000 # convert from kJ/(kg.K) to J/(kg.K)

#		fluid enthalpy
		h = IAPWS97(P = Pdata, T = Tdata).h
#		heat transfer coefficient
		x_tt = 1
	    
	result[0] = rho*result[0]
	result[1] = k*result[1]
	result[2] = mu*result[2]
	result[3] = cp*result[3]
	result[4] = h*result[4]
	result[5] = x_tt*result[5]
	
	return result

	
