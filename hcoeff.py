# module for calculate heat transfer coefficient 

# lib to use
import math
from scipy import optimize as so
import numpy as np
import fluid_data as fd
from iapws import IAPWS97

# constant defination
PI = 3.141592653
kw = 18.3 # wall conductivity in W/(m*K) 
#kw = 10.3 # wall conductivity in W/(m*K) 
PD_r = 1.3 #pitch diameter ration for boundle
g = -9.81 # gravity constant

# global statement 
# T for temperature in K
# p for Pressure in MPa
# x for steam gas fraction

# function for Reynold number calculation
def Reynold(dh,mdot,A,mu):  #eq diameter in m, mass flow in kg/s, flow area in m^2, dynamic visscosity in pa.s
	result = mdot*dh/(mu*A)

	return result

# function for Prandtl number calculation
def Prandtl(cp,mu,k): #specific heat capacity in J/(kg*K), dynamic visscosity in pa.s, heat conductivity in W/(m*K) 
	result = cp*mu/k

	return result

# function for Grashof number
def Grashof(Pdata,Tdata,Tw,dh):
	# difference between bulk and wall
	dT = abs(Tw - Tdata)
	# film temperature and density
	T_film = 1/2 * (Tw + Tdata)
	rho_film = IAPWS97(P = Pdata, T = T_film).rho
	# bulk thermal expansion and viscosity
	beta_l = IAPWS97(P = Pdata, T = Tdata).alfav
	mu_l = IAPWS97(P = Pdata, T = Tdata).mu

	result = abs(g)*beta_l*dT*dh**3/(mu_l/rho_film)**2

	return result



# function for colebrook
def f_colebrook(Re):
	def ff(Re): return lambda x: 1/math.sqrt(x) + 2.0*math.log10(2.51/(Re*math.sqrt(x))) #equition for solve f_colebrook
	
	if Re <= 2300:
		result = 64/Re
	elif (Re >= 2300) & (Re <= 1.0e4):
	        result = so.fsolve(ff(Re), 0.0306)
	elif (Re >= 1.0e4) & (Re <= 1.0e5):
	        result = so.fsolve(ff(Re), 0.024)
	elif (Re >= 1.0e5) & (Re <= 1.0e6):
	        result = so.fsolve(ff(Re),0.015)
	elif (Re >= 1.0e6) & (Re <= 1.0e7):
	        result = so.fsolve(ff(Re), 0.009)
	else:
	        result = so.fsolve(ff(Re), 0.005)
	
	return result

# function for Filonenko
def f_Filonenko(Re):
	result = (1.58*math.log(Re) - 3.28)**(-2)

	return result

# function for Nusselt number, single phse Gnielinski correction
def Nu_sp(Re,Pr):
#	f = f_colebrook(Re)
#	result = (f/8)*(Re-1000)*Pr / (1+12.7*(f/8)**(0.5) * (Pr*(2/3)-1))
	f = f_Filonenko(Re)
	result = (f/2)*(Re-1000)*Pr / (1+12.7*(f/2)**(0.5) * (Pr*(2/3)-1))
    
	return result

# function for Nusselt number, single phse Natrual Circulation
def Nu_NC_turb(Gr,Pr):
	
	result = 0.1*(Gr*Pr)**(1/3)

	return result

# function for Nusselt number with boundle correction
def Nu_Boundle(Re,Pr):
	C_eg = 0.028 * PD_r -0.006

	result = C_eg * Re**0.8 * Pr**0.33

	return result

# function for Nusselt number of subcooling with Saha and Zuber correlation
#def Nu_subcool(Re,Pr,dh,kl,qflx,h_sp):
def Nu_subcool(Pdata,Tdata,Tw):
	# dT between boiling Temperature and bulking Temperature
#	dT = 0.00135 * qflx/h_sp * Re**(1/2)
#	dT = qflx/(5*h_sp)
	T_sat = IAPWS97(P = Pdata, x = 0).T

#	result = qflx*dh/(kl*dT)
	result = 455*(T_sat - Tdata)/(Tw - Tdata)

	return result

# heat transfer coefficient for subcooled boiling
#def h_subcool(Pdata,Tdata,Tw,Tonb,h_sp,dh,kl):
def h_subcool(Pdata,Tdata,Tw,Tonb,h_sp):

	h0 = 5600 # in unit (W/m2*K)
	qflx0 = 20000 # in unit W/m^2

	# boling tempature 
	T_sat = IAPWS97(P = Pdata, x = 0).T
	# wall bulk temperature difference
	dT_bulk = Tw - Tdata
	# wall boling temperature difference
	dT_sat = Tw - T_sat
	# on set of nucleate boiling temperature difference
	dT_onb = Tonb - T_sat

	# single phase forced concection heat flux
	qflx_FC = h_sp * dT_bulk # in W/m^2

	# reduced pressure
	P_r = IAPWS97(P = Pdata, T = Tdata).Pr
	
	# correction factor 
	Fp = 1.73*P_r**0.27+(6.1+0.68/(1-P_r)) * P_r**2
	# power factor 
	n = 0.9 - 0.3 * P_r**0.15

	# pool boiling heat flux
	qflx_PB = (h0*Fp/qflx0**n)**(1/(1-n)) * dT_sat**(1/(1-n))

	# onb heat flux
	qflx_BI = (h0*Fp/qflx0**n)**(1/(1-n)) * dT_onb**(1/(1-n))	

	# nucleate boling heat flux
#	qflx_NB = (qflx_FC + (qflx_PB - qflx_BI)**3)**(1/3)
#	qflx_NB = qflx_FC + (qflx_PB**2 - qflx_BI**2)**(1/2)
	qflx_NB = qflx_FC + (qflx_PB**3 - qflx_BI**3)**(1/3)
#	qflx_NB = qflx_FC +(1/1+math.exp(-(Tw - Tonb)))* qflx_PB

#	print ('different Temperatures',Tdata,Tw,Tonb,T_sat)	
	print ('predict flux',qflx_NB,qflx_FC,qflx_PB,qflx_BI)	
	result = qflx_NB/dT_bulk

	print ('subcooled',result, h_sp)
	return result

# heat transfer coefficient in subcooled region 
def h_sp(Nu,k,dh):
	result = Nu*k/dh
    
	return result

# heat transfer coefficnt for pool boiling by cooper correlation
def h_pool(Pdata,xdata,qflx):
#	calculate reduced pressure
	P_r = IAPWS97(P = Pdata, x = xdata).Pr
#	define the Molecular Weight, in g/mol
	M = 18
#	assume the wall is smooth wall
	Rp = 1 # in um
#	power order of reduced pressure
	m_pr = 0.12 - math.log10(Rp)

#	calculate the pool boiling heat transfer coefficient
	result = 55 * P_r**m_pr * (-math.log10(P_r))**(-0.55) * M**(-0.5) * qflx**0.67

	return result

# heat transfer coefficient in saturate region, Dittus-Boelter correlation
def h_spDB(Re, Pr, k, dh):

	result = 0.023 * k/dh * Re**(0.8) * Pr**(0.4) 
    
	return result    

# heat transfer coefficient  in saturate region by winterton and gungor correlation
def h_satwg(pdata,xdata,mdot,Re,Pr,k,dh,x_tt,r,L,n_tube,qflx):

#	single phase heat transfer coefficient
	Nu = max(8.0,Nu_sp(Re,Pr))
#	h_sp = h_sp(Nu,k,dh)
	h_sp = h_spDB(Re,Pr,k,dh)
#	pool boiling heat transfer coefficient
	h_nc = h_pool(pdata,xdata,qflx)

#	latent heat 
	h_fg = IAPWS97(P = pdata, x=1).h - IAPWS97(P = pdata, x=0).h # in kJ/kg
#	calculate boiling number
	qflx = qflx/1000 # convert unit from W/m^(-2) to kW/m^(-2)
	Bo = qflx/(mdot*h_fg)

#	correlation factor
	E = 1 + 24000*Bo*1.16+1.37*(1/x_tt)**0.86
	S = 1/(1+1/1.15e-6*E**2*Re**1.17)

	result = E*h_sp + S*h_nc
#	print (result)	
	return result

# heat transfer coefficient  in saturate region by winterton and gungor correlation
def h_satwg2(pdata,xdata,mdot,Re,Pr,k,dh,A_xsc,r,L,n_tube,qflx):

	# steam and liquid density
	rho_g = IAPWS97(P = pdata, x = 1).rho
	rho_l = IAPWS97(P = pdata, x = 0).rho
	# single phase heat transfer coefficient
	Nu = max(8.0,Nu_sp(Re,Pr))
	h_l = h_sp(Nu,k,dh)
#	h_sp = h_spDB(Re,Pr,k,dh)
	
	# latant heat
	h_fg = IAPWS97(P = pdata, x=1).h - IAPWS97(P = pdata, x=0).h # in kJ/kg
	# calculate boiling number
	qflx = qflx/1000 # convert unit from W/m^(-2) to kW/m^(-2)
	Bo = qflx/(mdot*h_fg)
        # calculate Fr number
	G = mdot/A_xsc
	Fr = G**2/(g*dh*rho_l**2)

	# correlation factor
	S = 1+3000*Bo**0.86
	F = 1.12*(xdata/(1-xdata))**0.75*(rho_l/rho_g)**0.41	

	result = (F+S)*h_l

	return result


# heat transfer coefficient  in saturate region by Sun-Mishima correlation
def h_satSM(Pdata,xdata,mdot,k,dh,A_xsc,qflx):
	# steam and liquid density
	rho_g = IAPWS97(P = Pdata, x = 1).rho
	rho_l = IAPWS97(P = Pdata, x = 0).rho
	# steam and liquid viscosity
	mu_g = IAPWS97(P = Pdata, x = 1).mu
	mu_l = IAPWS97(P = Pdata, x = 0).mu
	# steam surface tension
	sigma = IAPWS97(P = Pdata, x = xdata).sigma
	# calculate ratio of  two phase viscosity and liquid viscosity , based on lin et la
	mu_tw = mu_g*mu_l/(mu_g + xdata**1.4*(mu_l - mu_g))
	# calculate liquid-only Reynold number
	Re_lo = Reynold(dh,mdot,A_xsc, mu_tw)
#	latent heat 
	h_fg = IAPWS97(P = Pdata, x=1).h - IAPWS97(P = Pdata, x=0).h # in kJ/kg
	# calculate boiling number
	qflx = qflx/1000 # convert unit from W/m^(-2) to kW/m^(-2)
	Bo = qflx/(mdot*h_fg)
	# calculate liquid Weber number
	mflx = mdot/A_xsc 
	We_l = mflx**2*dh/(rho_l*sigma)

	result = 6* Re_lo**1.05*Bo**0.54/(We_l**0.191*(rho_l/rho_g)**0.142) * k/dh

	return result

# heat transfer area
def A_transfer(r,L):
	result = 2*PI*r*L

	return result

# heat flux
def q_flx(q,r,L,n_tube):
	A_wall = A_transfer(r,L)

	q = q*1000/n_tube # convert the heat source from kW to W and divied by number of tubes
	result = q/A_wall # in unit W/m-2

	return result

# Temperature of onset of nucleate boiling
def T_onb(Tdata,Pdata,h_l,k_l):
	phi = 38 # contact angle for stainless steel, in degree
	phi = phi * PI/180 # convert degree to rradian
	# correction factor of contact angle
	F_phi = 1 - math.exp(-phi**3 - 0.5*phi)
	# boiling Tempeature
	T_sat = IAPWS97(P = Pdata, x = 0.0).T
	# serface tension
	sigma = IAPWS97(P = Pdata, T = Tdata).sigma
	# latent heat 
	h_fg = IAPWS97(P = Pdata, x=1).h - IAPWS97(P = Pdata, x=0).h # in kJ/kg
	h_fg = h_fg * 1000 # convert to J/kg
	# vapor density
	rho_g = IAPWS97(P = Pdata, x = 1.0).rho

	dT_onb_sat = 2*h_l/F_phi**2 * sigma*T_sat/(rho_g*h_fg*k_l)
	dT_sub = T_sat - Tdata

	# calculate the Temperature of onset of nucleate boiling
	result = Tdata+1/4 * (math.sqrt(dT_onb_sat)+math.sqrt(dT_onb_sat+4*dT_sub))**2

	return result

# Temperature of the secondary side wall
def Tw(mdot_pri,A_pri,r_pri,r_sec,L,dh,n_tube,T_pri,P_pri,q):
	# array for both sides temperature
	result = np.ones(2)
	# calcualate single liquid phase fluid property
	k_l = IAPWS97(P = P_pri, T = T_pri).k	
	cp_l = IAPWS97(P = P_pri, T = T_pri).cp * 1000 # convert from kJ/(kg K) to J/(kg K)
	mu_l = IAPWS97(P = P_pri, T = T_pri).mu	
	
	Re = max(100,Reynold(dh, mdot_pri, A_pri ,mu_l)) # calcualate single liquid phase Reynold number
	Pr = Prandtl(cp_l,mu_l,k_l) # calcualate single liquid phase Prandtl number

	F_bundle = 1 + 0.9120* Re**(-0.1)*Pr**(0.4)*(1-2.0043*math.exp(-PD_r)) # Factor for Bundle geometry, Finite array
	Nu = max(4.36,Nu_sp(Re, Pr)) * F_bundle  
	# primary heat transfer coefficient	
	h = h_sp(Nu, k_l, dh)

	# heat transfer area
	A_wall = 2*PI*r_pri * L

	# calculate the heat flux through the wall
	q = q*1000/n_tube # convert the heat source from kW to W and divied by number of tubes

	# primary side tube wall temperature
	Tw_pri = T_pri - q/(h*A_wall)

	Rw = (math.log(r_pri/r_sec)/(2*PI*L*kw)) # wall Resistance
	
	# secondary side tube wall temperature
	Tw_sec = Tw_pri - q*Rw

	result[0] = Tw_pri*result[0]
	result[1] = Tw_sec*result[1]

	return result
	
# heat transfer coefficient between fluid and wall
#def UA(mdot, A_xsc, r1, r2, L, dh, Tdata, Pdata,x_m, n_tube, q):
def UA(mdot, A_xsc, r1, r2, L, dh, Tdata, Pdata,x_ip1,x_i, n_tube, q, Tw):

# parameter statement:
#	mdot for mass flow rate in kg/s
#	r1 for heat transfer radius in m
#	r2 for another radius of the tube
#	L for cell length in m

#	average vapor quality
	x_m = 1/2*(x_ip1 + x_i)		
#	heat transfer area
	A_wall = 2*PI*r1 * L
#	calcualate two phase fluid property	
	fdata_tw = fd.fdata(Tdata,Pdata,x_m) 
	x_tt = fdata_tw[5]


#	calculate the heat flux through the wall
	q = q*1000/n_tube # convert the heat source from kW to W and divied by number of tubes
	q_flx = q/A_wall # in unit W/m-2
#	single phase heat transfer coefficeint 
	if (x_m == 0) or (x_m == 1.0):
	#	calcualate single liquid phase fluid property
		k_l = IAPWS97(P = Pdata, T = Tdata).k	
		cp_l = IAPWS97(P = Pdata, T = Tdata).cp * 1000 # convert from kJ/(kg K) to J/(kg K)
		mu_l = IAPWS97(P = Pdata, T = Tdata).mu	
	
		Re = max(100,Reynold(dh, mdot, A_xsc ,mu_l)) # calcualate single liquid phase Reynold number
		Pr = Prandtl(cp_l,mu_l,k_l) # calcualate single liquid phase Prandtl number

		if r1 >= r2: # primary side
		#	calcualate single liquid phase fluid property near the wall
			k_w = IAPWS97(P = Pdata, T = Tw).k	
			cp_w = IAPWS97(P = Pdata, T = Tw).cp * 1000 # convert from kJ/(kg K) to J/(kg K)
			mu_w = IAPWS97(P = Pdata, T = Tw).mu	

			Pr_w = Prandtl(cp_w,mu_w,k_w) # calcualate single liquid phase Prandtl number at wall
			
			F_correct = (Pr/Pr_w)**0.11
	#		F_bundle = 0.9217 + 0.1478*PD_r - 0.1130*math.exp(-7*PD_r - 1) # Factor for Bundle geometry, Infinate array
			
			F_bundle = 1 + 0.9120* Re**(-0.1)*Pr**(0.4)*(1-2.0043*math.exp(-PD_r)) # Factor for Bundle geometry, Finite array
	#		print (F_bundle,F_bundle_2)
			Nu = max(4.36,Nu_sp(Re, Pr)) * F_bundle# * F_correct
			
			h = h_sp(Nu, k_l, dh) 
		else: # secondary side
		#	calcualate single liquid phase fluid property near the wall
			k_w = IAPWS97(P = Pdata, T = Tw).k	
			cp_w = IAPWS97(P = Pdata, T = Tw).cp * 1000 # convert from kJ/(kg K) to J/(kg K)
			mu_w = IAPWS97(P = Pdata, T = Tw).mu	
		
			Pr_w = Prandtl(cp_w,mu_w,k_w) # calcualate single liquid phase Prandtl number at wall
		#	fluid property corrector
			F_correct = (Pr/Pr_w)**0.11
		
		#	natural circulation Nusult number
			Gr=  Grashof(Pdata,Tdata,Tw,dh)
			Nu_NC = Nu_NC_turb(Gr, Pr)
		#	force convection Nusult number	
			Nu_FC = max(4.36, Nu_sp(Re, Pr))	
		
			Nu = max(Nu_FC, Nu_NC) * F_correct 
			h = h_sp(Nu,k_l,dh) 
#			h_l = h_sp(Nu,k_l,dh) 
#			h_FC = h_sp(Nu_FC,k_l,dh)
#			if x_m == 0:		
#				Tonb =  T_onb(Tdata,Pdata,h_l,k_l)
#				if Tw > Tonb:
#					print ('ONB!')
#					print (Tonb, Tw)
#				#	h = h_subcool(Pdata,Tdata,Tw,Tonb,h_l) 
#				#	h = h_satSM(Pdata,x_m,mdot,k_l,dh,A_xsc,q_flx)
#				#	h = h_subcool(Pdata,479.2,528,524.23,h_l) 
#				#	h = h_sp(Nu, k_l, dh)
#					print ('heat transfer coeffcient',h)
#				else:
#					print ('Single!')
#				#	print ('Wall Temperature', Tw)
#				#	h = h_sp(Nu, k_l, dh) 
#			elif x_m == 1:
#				h = h_sp(Nu, k_l, dh)
#	two phase heat transfer coefficient
	elif x_m > 0  and x_m < 1.0:	
	#	calcualate single liquid phase fluid property
		k_l = IAPWS97(P = Pdata, x = 0.0).k	
		cp_l = IAPWS97(P = Pdata, x = 0.0).cp * 1000 # convert from kJ/(kg K) to J/(kg K)
		mu_l = IAPWS97(P = Pdata, x = 0.0).mu	

		Re = max(100,Reynold(dh, mdot, A_xsc ,mu_l)) # calcualate single liquid phase Reynold number
		Pr = Prandtl(cp_l,mu_l,k_l) # calcualate single liquid phase Prandtl number

#		h = h_satwg(Pdata,x_m, mdot, Re,Pr,k_l,dh,x_tt,r1,L,n_tube,q_flx)
#		h = h_satwg2(Pdata,x_m, mdot, Re,Pr,k_l,dh,A_xsc,r1,L,n_tube,q_flx)
		h = h_satSM(Pdata,x_m,mdot,k_l,dh,A_xsc,q_flx)
	#	print ('two phase heat transfer',h)

	Rw = (math.log(max(r1,r2)/min(r1,r2))/(2*PI*L*kw))/2 # half of the total wall conductance of cylinder
#	Rw = abs(r1-r2)/kw/2 # half of the total wall conductance of slab 
	result = 1/(Rw + 1/(h*A_wall))
    
	return result

# calculation of critical heat flux, biasi corrolation
def q_chf(mdot,Pdata,x_ip1,x_i,A_xsc,dh, h_m, h_boil):

#	boil length
	LB = h_m - h_boil
#	mass flux
	G = mdot/A_xsc
#	convert Pressure from MPa to Bar
	P = Pdata*10
#	average vapor quality
	x_m = 1/2 * (x_ip1+x_i)	

#	latent heat 
	h_fg = IAPWS97(P = Pdata, x=1).h - IAPWS97(P = Pdata, x=0).h # in kJ/kg
	h_fg = 1000*h_fg # convert to J/kg

#	F and H function
	F_P = 0.7249 + 0.099*P*math.exp(-0.032*P)
	H_P = -1.159+0.149*P*math.exp(-0.019*P) + 8.99*P/(10+P**2)

#	x1 parameters
	A1 = 1.0
	B1 = (1.048e-8 * G**(1.6) * dh**(1.4) * h_fg)/H_P
#	x2 parameters
	A2 = F_P/G**(1/6)
	B2 = 5.707e-8 * G**(7/6) * dh**(1/4) * h_fg
#	x values
	x1 = A1*LB/(B1+LB)
	x2 = A2*LB/(B2+LB)
#	print (x1,x2)
		
#	biasi corrolation
#	result = max(x1,x2)
	result = 15.048e7* (100*dh)**0.4 * G**(-0.6) * H_P * (1-x_m)	

	return result	
