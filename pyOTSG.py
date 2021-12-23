# Dynamic model of once-througn steam generator 

#imported library
import os.path
import sys
import numpy as np
from iapws import IAPWS97
import math
import hcoeff as hc
import pdrop as pd
import fluid_data as fd
import eqsolver as eq
import time
 
# steady state output file name
steadyfile = 'steady_output.txt'

# transient state output file name
transfile = 'trans_output.txt'

# delete transient result
if os.path.isfile(transfile):
	os.remove(transfile)
	os.remove('Transfered_heat.txt')
	os.remove('InoutT.txt')
	os.remove('InoutP.txt')
	os.remove('Phasechange_Loc.txt')
	os.remove('udiff.txt')
if os.path.isfile('Massflow_compare.txt'):
	os.remove('Massflow_compare.txt')
if os.path.isfile('enthalpy_compare.txt'):
	os.remove('enthalpy_compare.txt')
if os.path.isfile('iteration.txt'):
	os.remove('iteration.txt')

# simulation mode (1 for run steady and transient, 0 for run jump steady, 2 for steady only)
mode = 1

# define convergence 
eps = 1e-4
eps_trans = 1e-4

# define number of nodes
n = 311 
# define number of loops
n_loop = 100

# define lib boundary
#dp_lim = 10
# define phase change boundary
x_b = 1.0

#################### transient information #####################
# start and end time 
time_start = 0
time_end = 120

# define timestep
dt = 0.05 # in s

# steam generator reaction time, in s
t_trip = 5
t_iso = 5  
###############################################################

# define constant
T_ck = 273.15 # convert factor from degree C to degree K
PI = 3.141592653
#kw = 18.3 #W/(m.K), mean conductivity of Inconel 600
kw = 10.3 #W/(m.K), mean conductivity of Inconel 600

# defination of SG geometry
l_tube = 610 # in cm, the length of SG 
#l_tube = 550 # in cm, the length of SG 

d_tube = 0.12 # in cm, thickness of  tube
od_tube = 1.96+2*d_tube # in cm, outer diameter of secondary tube

CSA_SG = 2.751 #in m^2, cross section area of SG
PD_r = 1.3 # PD_ratio

#defination of SG parameters
mf_pri = 995.36 # in kg/s, mass flow rate of primary circuit
T_inlet_pri = 309  # in degree C, inlet temperature of primary circuit
P_pri = 14.5 # in MPa, system pressure of primary circuit 


mf_sec = 94 # in kg/s, mass flow rate of secondary circuit
#mf_sec = 83.15 # in kg/s, mass flow rate of secondary circuit
#mf_sec = 86.8 # in kg/s, mass flow rate of secondary circuit
T_fw = 180  # in degree C, inlet temperature of secondary circuit
#T_fw = 138  # in degree C, inlet temperature of secondary circuit
#T_fw = 120  # in degree C, inlet temperature of secondary circuit
T_steam = 285  # in degree C, outlet temperature of secondary circuit
P_sec = 4.0 # in Mpa, system pressure of secondary circuit 


# convert degree C to degree K
T_inlet_pri = T_inlet_pri + T_ck
T_fw = T_fw + T_ck
T_steam = T_steam + T_ck

#convert data in cm to m
l_tube = l_tube/100
d_tube = d_tube/100
od_tube = od_tube/100

# primary and secondary flow area calculation
#primary side
r_pri = od_tube/2 # in m, tube radius
A_pri = ((PD_r**2 - PI/4)/PD_r**2) * CSA_SG # in m^2, primary area

#secondary side
r_sec = od_tube/2 - d_tube # in m, tube radius
A_sec = ((PI/4)/PD_r**2) * CSA_SG  # in m^2, secondary area

#number of tubes
n_tube = int(float(A_sec/(r_sec**2*PI))) # number of tubes of SG

#equivalent diameter calculation
dh_pri = 4*((PD_r * od_tube)**2 - PI/4 *od_tube**2)/(PI*od_tube) #in m, primary side equivalent diameter
dh_sec = od_tube - 2*d_tube  #in m, secondary side equivalent diamter

#cell total length 
l_cell = l_tube/(n-1) 
# initial fluid data calculation 

#enthalpy of primary coolant inlet
h_inlet_pri = IAPWS97(P = P_pri, T = T_inlet_pri).h # in kJ/kg, same below

# enthalpy  of secondary feedwater to superheated steam, in kJ/kg
h_fw_sec = IAPWS97(P = P_sec,T = T_fw).h #feedwater inlet enthalpy
h_boil_sec = IAPWS97(P = P_sec, x = 0).h #boiling point enthalpy
h_sat_sec = IAPWS97(P = P_sec, x = 1).h #satured steam enthalpy
h_steam_sec = IAPWS97(P = P_sec, T = T_steam).h #outlet steam enthalpy

h_latent_sec = h_sat_sec - h_boil_sec
# secondary side two phase flow temperature, in K
T_tp_sec = IAPWS97(P = P_sec, x = 1).T

# mass flow rate of primary and secondary side, in kg/s
G_pri = mf_pri
G_sec = mf_sec

#========================= basic calculation for boundary =================================

#calculation of total power
h_s1 = h_fw_sec
h_s3 = h_boil_sec
h_s5 = h_sat_sec
h_s7 = h_steam_sec

Q_s2 = G_sec*(h_s3-h_s1)  # subcooling region, in kW
Q_s4 = G_sec*(h_s5-h_s3)  # boiling region, in kW
Q_s6 = G_sec*(h_s7-h_s5)  # superheated region, in kW

# total power 
Q = Q_s2 + Q_s4 + Q_s6 #in kW
#print ('total power',Q/1000,'[MW]')

# initial transfered heat array
q_rough = np.ones(n-1)

# primary outlet temperature 
h_outlet_pri = h_inlet_pri - Q/G_pri # primary side outlet enthalpy 
T_outlet_pri = IAPWS97(P = P_pri,h = h_outlet_pri).T
T_outlet_pri = 271 + T_ck

# estimate of primary pressure drop
T_ave = 1/2 * (T_outlet_pri + T_inlet_pri)
rho_pri_in = IAPWS97(P = P_pri,h = h_inlet_pri).rho
rho_pri_out = IAPWS97(P = P_pri,h = h_outlet_pri).rho
Dp_pri_est =  pd.pd_single(G_pri,G_pri,T_ave,P_pri,rho_pri_in,rho_pri_out,0.0,dh_pri,A_pri,l_cell,-1)
#========================= check steady state result =================================
if os.path.isfile(steadyfile):
	print ('steady result exist')
	if mode == 1:
		mode = 0
	if mode == 2:
		os.remove(steadyfile)
		print ('removed')
else:
	print ('no steady state results found')
	mode = 1

if mode == 1 or mode == 2:
	print ('\033[1m'+'**** start stead state calculation ****'+ '\033[0m')
#======================== steady state calculation =================================
#	initial temperature array
	Tarray_pri = T_outlet_pri*np.ones(n)
#	Tarray_pri = T_inlet_pri*np.ones(n)
	Tarray_sec = T_fw*np.ones(n)

#	initial wall 
#	Tarray_wall = np.ones(n)
#	Tarray_wall = ((T_outlet_pri + T_fw)/2) * Tarray_wall

#	initial pressure array
	Parray_pri = (P_pri - Dp_pri_est)*np.ones(n)
	Parray_sec = P_sec*np.ones(n)

#	initial pressure array
	Garray_pri = G_pri*np.ones(n)
	Garray_sec = G_sec*np.ones(n)

#	initial density array
	rho_pri = IAPWS97( P = P_pri,T = T_outlet_pri).rho
#	rho_pri = IAPWS97( P = P_pri,T = T_inlet_pri).rho
	rho_sec = IAPWS97( P = P_sec,T = T_fw).rho
	rhoarray_pri = rho_pri * np.ones(n)
	rhoarray_sec = rho_sec * np.ones(n)

#	initial velocity array
	uarray_pri = G_pri/(rhoarray_pri*A_pri)
	uarray_sec = G_sec/(rhoarray_sec*A_sec)
	
#	initial enthalpy array
	h_pri = h_outlet_pri*np.ones(n)
#	h_pri = h_inlet_pri*np.ones(n)
	h_sec = h_fw_sec*np.ones(n)
	
#	initial vaper fraction array
	x_pri = np.zeros(n)
	x_f = np.zeros(n)

#########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#############
	dp_pri_array = []
	dp_sec_array = []	
#########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#############
#	initalize boiling condition, 0 for liquid, 1 for boiling, 2 for superheated
	sat = 0
	print ('**********************************')	
	print ('** iteration for rough estimate **')	
	print ('**********************************')	
	
	print ('***single phase liquid region***')
	print('node Temperature', 0 ,Tarray_pri[0],Tarray_sec[0], x_f[0])
	for i in range(n-1):
#		initial guess of total heat in each region
		if sat == 0:
			q = Q_s2/((n-1)/4)
		elif sat == 1:
			q = Q_s4/((n-1)/3)
			#q = Q_s4/40
		else:
			q = Q_s6/((n-1)/3)
	
#		guess of  enthalpy of each side
		h_pri[i+1] = h_pri[i] + q/G_pri	
		h_sec[i+1] = h_sec[i] + q/G_sec

#		guess of next node Temperature
		Tarray_pri[i+1] = Tarray_pri[i] + (T_inlet_pri-T_outlet_pri)/n
		if sat == 0 or sat == 2:
			Tarray_sec[i+1] = Tarray_sec[i] + (T_steam - T_fw)/n
		else:
			Tarray_sec[i+1] = IAPWS97(P = Parray_sec[i+1], h = h_sec[i+1]).T

#		calculate cell mean temperature of each side
		T_pri_m = 1/2 * (Tarray_pri[i] + Tarray_pri[i+1])
		T_sec_m = 1/2 * (Tarray_sec[i] + Tarray_sec[i+1])

#		guess of next node Pressure
		Parray_pri[i+1] = Parray_pri[i] + pd.pd_single(G_pri,G_pri,T_pri_m,Parray_pri[i],rhoarray_pri[i+1],rhoarray_pri[i],0,dh_pri,A_pri,l_cell,1)
		Parray_sec[i+1] = Parray_sec[i] - pd.pd_single(G_sec,G_sec,T_sec_m,Parray_sec[i],rhoarray_sec[i+1],rhoarray_sec[i],0,dh_sec,A_sec,l_cell,1)

#		guess of next node density
		rhoarray_pri[i+1] = IAPWS97(P = Parray_pri[i+1], h = h_pri[i+1]).rho
		rhoarray_sec[i+1] = IAPWS97(P = Parray_sec[i+1], h = h_sec[i+1]).rho

#		guess of next node velocity
		uarray_pri[i+1] = G_pri/(rhoarray_pri[i+1]*A_pri)
		uarray_sec[i+1] = G_sec/(rhoarray_sec[i+1]*A_sec)	

#		guess of next node vapor fraction
		x_f[i+1] = IAPWS97(P = Parray_sec[i+1], h = h_sec[i+1]).x

#		guess of converge data
		h_pri_old = h_pri[i+1] 
		h_sec_old = h_sec[i+1] 
		u_pri_old = uarray_pri[i+1] 
		u_sec_old = uarray_sec[i+1] 
		rho_pri_old = rhoarray_pri[i+1] 
		rho_sec_old = rhoarray_sec[i+1]
		q_old = q 

#		node data iteration	
		for j in range(n_loop):
	#		cell average temperature
			T_pri_m = 1/2 * (Tarray_pri[i] + Tarray_pri[i+1])
			T_sec_m = 1/2 * (Tarray_sec[i] + Tarray_sec[i+1])
	#		cell average pressure
			P_pri_m = 1/2 * (Parray_pri[i] + Parray_pri[i+1])
			P_sec_m = 1/2 * (Parray_sec[i] + Parray_sec[i+1])
	#		cell average mass flow
			G_pri_m = 1/2 * (Garray_pri[i] + Garray_pri[i+1])
			G_sec_m = 1/2 * (Garray_sec[i] + Garray_sec[i+1])
	#		cell average vapor fraction
			x_pri_m = 0
			x_sec_m = 1/2 * (x_f[i] + x_f[i+1])
	
	#		slove momentum equation, mass continuty and energy equation
	#		solve mass continuity equation
	#		momentum equation
	#		total pressure drop by corrolations
			Parray_pri[i+1] = Parray_pri[i] + pd.pd_single(Garray_pri[i+1],Garray_pri[i],T_pri_m,P_pri_m,rhoarray_pri[i+1],rhoarray_pri[i],0,dh_pri,A_pri,l_cell,1)
	
			if x_sec_m == 0 or x_sec_m >= x_b:
	#			single phase pressure drop
				Parray_sec[i+1] = Parray_sec[i] - pd.pd_single(Garray_sec[i+1],Garray_sec[i],T_sec_m,P_sec_m,rhoarray_sec[i+1],rhoarray_sec[i],x_sec_m,dh_sec,A_sec,l_cell,1)
			elif x_sec_m > 0 and x_sec_m < x_b:
				Parray_sec[i+1] = Parray_sec[i] - pd.pd_twHEM(Garray_sec[i+1],Garray_sec[i],T_sec_m,P_sec_m,rhoarray_sec[i+1],rhoarray_sec[i],x_f[i+1],x_f[i],dh_sec,A_sec,l_cell,1)
			
	#		solve energy equation
	#		calculate heat source or sink for the cell
	#		log mean temperature
			dT_up = Tarray_pri[i+1] - Tarray_sec[i+1]
			dT_low = Tarray_pri[i] - Tarray_sec[i]
			dT_m = (dT_up - dT_low)/math.log(dT_up/dT_low)
	#		calcualte secondary wall temperature
			Twall = hc.Tw(G_pri_m,A_pri,r_pri,r_sec,l_cell,dh_pri,n_tube,T_pri_m,P_pri_m,q)
	#		print (Twall)
	#		calculate overall heat transfer coefficient
			UA_pri = hc.UA(G_pri_m,A_pri,r_pri,r_sec,l_cell,dh_pri,T_pri_m,P_pri_m,0.0,0.0, n_tube, q, Twall[0])
			UA_sec = hc.UA(G_sec_m,A_sec,r_sec,r_pri,l_cell,dh_sec,T_sec_m,P_sec_m,x_f[i+1],x_f[i], n_tube, q, Twall[1])
			UA = 1/(1/UA_pri + 1/UA_sec)
	#		print (UA)
	#		waring for nonsence value	
			if UA <= 0:
				print ('************* WARNING **************')
				print ('**negative heat transfer coeffient**')
				print ('************************************')
	#		calculate transfered heat
			q = float('%.10f'%(UA*dT_m/1000 * n_tube)) # in unit kW
		#	print ('transfered heat',q/1000)
		#	print ('\n')
		#	if j >= 30:
		#	if i >= 8:
		#		print (i)
		#		sys.exit()
	#		enthalpy primary side
			h_pri[i+1] = float('%.10f'%(h_pri[i] + q/G_pri_m))
	#		enthalpy secondary side
			h_sec[i+1] = float('%.10f'%(h_sec[i] + q/G_sec_m))
	#		recalculate the Temperature
			Tarray_pri[i+1] = float('%.10f'%IAPWS97( P = Parray_pri[i+1], h = h_pri[i+1]).T)
			Tarray_sec[i+1] = float('%.10f'%IAPWS97( P = Parray_sec[i+1], h = h_sec[i+1]).T)
			
	#		calculate the vapor quality
			x_f[i+1] = float('%.10f'%(IAPWS97( P = Parray_sec[i+1], h = h_sec[i+1]).x))

	#		calculate the density
			rhoarray_pri[i+1] = float('%.10f'%IAPWS97(P = Parray_pri[i+1], h = h_pri[i+1]).rho)
			rhoarray_sec[i+1] = float('%.10f'%IAPWS97(P = Parray_sec[i+1], h = h_sec[i+1]).rho)

	#		solve mass contiunity equation for new volocity
			uarray_pri[i+1] = float('%.10f'%(G_pri/(rhoarray_pri[i+1]*A_pri)))
			uarray_sec[i+1] = float('%.10f'%(G_sec/(rhoarray_sec[i+1]*A_sec)))

			Garray_pri[i+1] = float('%.10f'%(rhoarray_pri[i+1] * uarray_pri[i+1] * A_pri))
			Garray_sec[i+1] = float('%.10f'%(rhoarray_sec[i+1] * uarray_sec[i+1] * A_sec))

	#		recalculate the boiling enthalpy
			h_boil = float('%.10f'%(IAPWS97(P = Parray_sec[i+1], x = 0).h))
	#		recalculate the saturated enthalpy
			h_sat = float('%.10f'%(IAPWS97(P = Parray_sec[i+1], x = 1).h))
			if sat == 1:
	#			recalculate the critical quality
				h_m = (i+1/2) * l_cell
				q_chf = hc.q_chf(G_sec,P_sec_m,x_f[i],x_f[i+1],A_sec,dh_sec,h_m,height_boil)
				q_flx = hc.q_flx(q,r_sec,l_cell,n_tube)
				if q_flx > q_chf:
					print ('chf',i+1)
					print (q_chf,q_flx)
				#	sys.exit()
	
		#	if j >= n_loop*4/5:
			#	print (h_sec[i+1], h_sec_old)
		#		print (h_pri[i+1] - h_pri_old,h_sec[i+1] - h_sec_old)
		#		print (rhoarray_sec[i+1] - rho_sec_old)
		#		print (uarray_sec[i+1] - u_sec_old)
		#		print (q_old - q)

	#		check convergency			
			if (abs(rho_pri_old - rhoarray_pri[i+1]) < eps*10 and abs(rho_sec_old - rhoarray_sec[i+1]) < eps*10
				and abs(u_pri_old - uarray_pri[i+1]) < eps and abs(u_sec_old - uarray_sec[i+1]) < eps
					and abs(h_pri_old - h_pri[i+1]) < eps*10 and abs(h_sec_old - h_sec[i+1]) < eps*10
					    and abs(q_old - q) < eps*100):
				q_rough[i] = q
				break
			elif j>= (n_loop-1):
				print ('residual value', )
				print (abs(rho_pri_old - rhoarray_pri[i+1]) - eps*1e3,  abs(rho_sec_old - rhoarray_sec[i+1])- eps*1e3)
				print (abs(u_pri_old - uarray_pri[i+1]) - eps , abs(u_sec_old - uarray_sec[i+1]) - eps)
				print (abs(h_pri_old - h_pri[i+1]) - eps*10 , abs(h_sec_old - h_sec[i+1]) - eps*10)
				print (x_f[i],x_f[i+1])
				print (abs(q_old - q) - eps)
				sys.exit('unconverge')
	
	#		store current iteration value
			rho_pri_old = rhoarray_pri[i+1]
			rho_sec_old = rhoarray_sec[i+1]
		
			u_pri_old = uarray_pri[i+1]
			u_sec_old = uarray_sec[i+1]
	
			h_pri_old = h_pri[i+1]
			h_sec_old = h_sec[i+1]

			q_old = q
	
	#	determine whether boiling occurs
		if x_f[i+1] > 0 and x_f[i] == 0:
			sat = 1
	#		node of boiling and height of boilling
			height_boil = (i+1)*l_cell - l_cell/2

			print ('*****************')
			print ('boiling')
			print ('*****************')
			print ('node boiling', i+1, 'at height', '%.10f'%height_boil, '[m]')
	
	#	determine saturate steam
		elif x_f[i] < x_b and x_f[i+1] >= x_b:
	#		change liquid state 			
			sat = 2
			height_sat = (i+1)*l_cell - l_cell/2

			print ('*****************')
			print ('saturated')
			print ('*****************')
			print ('node saturated', i+1, 'at height','%.10f'%height_sat, '[m]')
	
#		output result		
		print ('node data', i+1,Parray_pri[i+1] ,Tarray_pri[i+1], Parray_sec[i+1], Tarray_sec[i+1],x_f[i+1])
		print (rhoarray_sec[i+1], uarray_sec[i+1],Garray_sec[i+1])	
#	total heat transfered
	Qtotal = sum(q_rough)
	Q_pri = (h_pri[n-1] - h_pri[0])*(G_pri/1000)
	Q_sec = (h_sec[n-1] - h_sec[0])*(G_sec/1000)
	
	print ('*****************************')	
	print (Qtotal/1000,Q_pri,Q_sec)
	print ('*****************************')	
			

	print ('*****************************')	
	print ('** counter flow correction **')	
	print ('*****************************')	
#	initial heat transfer
	q_trans = q_rough * np.ones(n-1)
	q_temp = q_rough * np.ones(n-1)
	Qtotal_old = sum(q_rough)
#	boundary condition correction
	Tarray_pri[n-1] = T_inlet_pri
	Parray_pri[n-1] = P_pri
	h_pri[n-1] = IAPWS97(P = P_pri, T = T_inlet_pri).h # in kJ/kg, same below

	for j in range(n_loop):	
	
	#	data primary side
		for i in range(n-1):
			u_pri_old = uarray_pri[n-i-1-1] 
			rho_pri_old = rhoarray_pri[n-i-1-1] 
			h_pri_old = h_pri[n-i-1-1]
			for k in range(n_loop):
			#	cell average data
				T_pri_m = 1/2 * (Tarray_pri[n-i-1] + Tarray_pri[n-i-1-1])
				P_pri_m = 1/2 * (Parray_pri[n-i-1] + Parray_pri[n-i-1-1])
				G_pri_m = 1/2 * (Garray_pri[n-i-1] + Garray_pri[n-i-1-1])
				x_pri_m = 1/2 * (x_pri[n-i-1] + x_pri[n-i-1-1])
		#		total pressure drop by corrolations
				Parray_pri[n-i-1-1] = Parray_pri[n-i-1] - pd.pd_single(Garray_pri[n-i-1-1],Garray_pri[n-i-1],T_pri_m,P_pri_m,rhoarray_pri[n-i-1-1],rhoarray_pri[n-i-1],x_pri_m,dh_pri,A_pri,l_cell,-1)
		
			#	slove  energy equation
				h_pri[n-i-1-1] = float('%.10f'%(h_pri[n-i-1] - q_trans[n-i-1-1]/G_pri_m))
			#	lookup data from IAPWS97
				Tarray_pri[n-i-1-1] = float('%.10f'%IAPWS97( P = Parray_pri[n-i-1-1], h = h_pri[n-i-1-1]).T)
				rhoarray_pri[n-i-1-1] = float('%.10f'%IAPWS97(P = Parray_pri[n-i-1-1], h = h_pri[n-i-1-1]).rho)
			#	solve mass continuty equation for velocity
				uarray_pri[n-i-1-1] = float('%.10f'%(G_pri_m/(rhoarray_pri[n-i-1-1]*A_pri)))
			#	recalculate mass flow rate
				Garray_pri[n-i-1-1] = float('%.10f'%(rhoarray_pri[n-i-1-1] * uarray_pri[n-i-1-1] * A_pri))
			
			#	check convergency for primary side			
				if (abs(rho_pri_old - rhoarray_pri[n-i-1-1]) < eps*10 and abs(u_pri_old - uarray_pri[n-i-1-1]) < eps):
		#				and abs(h_pri_old - h_pri[i+1]) < eps*10):
					break
				elif k >= (n_loop-1):
					print ('primary residual value',i+1 )
					print (abs(rho_pri_old - rhoarray_pri[i+1]) - eps*10)
					print (abs(u_pri_old - uarray_pri[i+1]) - eps)
					print (abs(h_pri_old - h_pri[i+1]) - eps*10)
					sys.exit('unconverge')
	
			#	store current iteration value
				rho_pri_old = rhoarray_pri[n-i-1-1]
				u_pri_old = uarray_pri[n-i-1-1]
				h_pri_old = h_pri[n-i-1-1]

	#	data secondary side
		for i in range(n-1):
			u_sec_old = uarray_sec[i+1] 
			rho_sec_old = rhoarray_sec[i+1]
			h_sec_old = h_sec[i+1]
			for k in range(n_loop):
			#	cell average data of primary side
				T_pri_m = 1/2 * (Tarray_pri[i] + Tarray_pri[i+1])
				P_pri_m = 1/2 * (Parray_pri[i] + Parray_pri[i+1])
				G_pri_m = 1/2 * (Garray_pri[i] + Garray_pri[i+1])
				x_pri_m = 1/2 * (x_pri[i] + x_pri[i+1])
			#	cell average data of secondary side
				T_sec_m = 1/2 * (Tarray_sec[i] + Tarray_sec[i+1])
				P_sec_m = 1/2 * (Parray_sec[i] + Parray_sec[i+1])
				G_sec_m = 1/2 * (Garray_sec[i] + Garray_sec[i+1])
				x_sec_m = 1/2 * (x_f[i] + x_f[i+1])
			
			#	totol presure drop by corrolation
				if x_sec_m == 0 or x_sec_m >= x_b:
				#	single phase pressure drop
					Parray_sec[i+1] = Parray_sec[i] - pd.pd_single(Garray_sec[i+1],Garray_sec[i],T_sec_m,P_sec_m,rhoarray_sec[i+1],rhoarray_sec[i],x_sec_m,dh_sec,A_sec,l_cell,1)
				elif x_sec_m > 0 and x_sec_m < x_b:
				#	two phase pressure drop
					Parray_sec[i+1] = Parray_sec[i] - pd.pd_twHEM(Garray_sec[i+1],Garray_sec[i],T_sec_m,P_sec_m,rhoarray_sec[i+1],rhoarray_sec[i],x_f[i+1],x_f[i],dh_sec,A_sec,l_cell,1)
	
			#	calculate heat source or sink for the cell
			#	log mean temperature
				dT_up = Tarray_pri[i+1] - Tarray_sec[i+1]
				dT_low = Tarray_pri[i] - Tarray_sec[i]
				dT_m = (dT_up - dT_low)/math.log(dT_up/dT_low)
			#	calcualte secondary wall temperature
				Twall = hc.Tw(G_pri_m,A_pri,r_pri,r_sec,l_cell,dh_pri,n_tube,T_pri_m,P_pri_m,q)
			#	calculate overall heat transfer coefficient
				UA_pri = hc.UA(G_pri_m,A_pri,r_pri,r_sec,l_cell,dh_pri,T_pri_m,P_pri_m,x_pri[i+1],x_pri[i], n_tube, q_trans[i],Twall[0])
				UA_sec = hc.UA(G_sec_m,A_sec,r_sec,r_pri,l_cell,dh_sec,T_sec_m,P_sec_m,x_f[i+1],x_f[i], n_tube, q_trans[i],Twall[1])
				UA = 1/(1/UA_pri + 1/UA_sec)
			#	waring for nonsence value	
				if UA <= 0:
					print ('************* WARNING **************')
					print ('**negative heat transfer coeffient**')
					print ('************************************')
			
			#	calculate transfered heat
				q = float('%.10f'%(UA*dT_m/1000 * n_tube)) # in unit kW
				q_trans[i] = q
			
			#	slove energy equation
				h_sec[i+1] = float('%.10f'%(h_sec[i] + q_trans[i]/G_sec_m))
			#	lookup data from IAPWS97
				Tarray_sec[i+1] = float('%.10f'%IAPWS97( P = Parray_sec[i+1], h = h_sec[i+1]).T)
				rhoarray_sec[i+1] = float('%.10f'%IAPWS97(P = Parray_sec[i+1], h = h_sec[i+1]).rho)
				x_f[i+1] = float('%.10f'%(IAPWS97( P = Parray_sec[i+1], h = h_sec[i+1]).x))
			
			#	solve mass countinuity equation
				uarray_sec[i+1] = float('%.10f'%(G_sec/(rhoarray_sec[i+1]*A_sec)))
			#	recalculate mass flow rate
				Garray_sec[i+1] = float('%.10f'%(rhoarray_sec[i+1] * uarray_sec[i+1] * A_sec))
	
			#	recalculate the boiling enthalpy
				h_boil = float('%.10f'%(IAPWS97(P = Parray_sec[i+1], x = 0).h))
			#	recalculate the saturated enthalpy
				h_sat = float('%.10f'%(IAPWS97(P = Parray_sec[i+1], x = 1).h))
			
			#	check convergency of secondary side	
				if (abs(rho_sec_old - rhoarray_sec[i+1]) < eps*10 and abs(u_sec_old - uarray_sec[i+1]) < eps
						and abs(h_sec_old - h_sec[i+1]) < eps*10 and abs(q_old - q) < eps*100):
					q_trans[i] = q
					break
				elif k >= (n_loop-1):
					print ('secondary residual value',i+1 )
					print (abs(rho_sec_old - rhoarray_sec[i+1])- eps*10)
					print (abs(u_sec_old - uarray_sec[i+1]) - eps)
					print (abs(h_sec_old - h_sec[i+1]) - eps*10)
					print (x_f[i],x_f[i+1])
					print (abs(q_old - q) - eps*10)
					sys.exit('unconverge')
		
			#	store current iteration value
				rho_sec_old = rhoarray_sec[i+1]
				u_sec_old = uarray_sec[i+1]
				h_sec_old = h_sec[i+1]
				q_old = q

	#	total heat transfered
		Qtotal = sum(q_trans)
		Q_pri = (h_pri[n-1] - h_pri[0])*(G_pri/1000)
		Q_sec = (h_sec[n-1] - h_sec[0])*(G_sec/1000)
		
		print (sum(q_trans)/1000,Q_pri,Q_sec)
			
	#	check convergency 
		con_check = abs(q_trans - q_temp)
		if all (item  < eps for item in con_check):
			print ('number of iterations', j)
			break	
		elif j >= (n_loop-1):
			print ('Heat residual value', )
			print (abs(Qtotal_old - Qtotal) - eps*10)
			sys.exit('unconverge')
		Qtotal_old = sum(q_trans)
		q_temp = q_trans * np.ones(n-1)

#	sys.exit()
#	redietermine the boiling and saturated point
	for i in range(n-1):
	#	determine boiling
		if x_f[i+1] > 0 and x_f[i] == 0:
	#		node of boiling and height of boilling
			height_boil = (i+1)*l_cell - l_cell/2
			print ('node boiling', i+1, 'at height', '%.10f'%height_boil, '[m]')
	
	#	determine saturate steam
		elif x_f[i] < x_b and x_f[i+1] >= x_b:
			height_sat = (i+1)*l_cell - l_cell/2
			print ('node saturated', i+1, 'at height','%.10f'%height_sat, '[m]')
#	rename array
	u_pri = uarray_pri * np.ones(n)
	u_sec = uarray_sec * np.ones(n)
	rho_pri = rhoarray_pri * np.ones(n)
	rho_sec = rhoarray_sec * np.ones(n)
	Q_transfer = q_trans * np.ones(len(q_trans))
	node = list(range(0,n))
	node = np.array(node)
	height = l_cell * node

#=========================steady state output file creation=================================
#	create node array
	node = list(range(0,n))
	node = np.array(node)
	height = l_cell * node
	data_pri = np.column_stack((height,Tarray_pri,h_pri))
	data_sec = np.column_stack((height,Tarray_sec,h_sec,x_f))
	cell_centre = height[0:n-1] + l_cell/2
	Q_trans = np.column_stack((cell_centre,q_trans))
	Parray_data = np.column_stack((height,Parray_pri,Parray_sec))
	u_data = np.column_stack((height,uarray_pri,uarray_sec))
	rho_data = np.column_stack((height,rhoarray_pri,rhoarray_sec))
	
	with open(steadyfile, 'w') as f:
		f.write('#'+'boiling point [m]    saturation point [m]\n')
		f.write(str('%.10f'%height_boil)+'    '+str('%.10f'%height_sat)+'\n')

		f.write('#'+'primary side data:\n')
		f.write('#'+'node height [m]    temperature [K]   enthalpy [kJ/kg] \n')
		for row in data_pri:	
			f.write('	'.join(map(str,row)) + '\n')
		f.write('\n')

		f.write('#'+'secondary side data:\n')
		f.write('#'+'node height [m]    temperature [K]   enthalpy [kJ/kg]   vapor fraction \n')
		for row in data_sec:	
			f.write('	'.join(map(str,row)) + '\n')
		f.write('\n')

		f.write('#'+'heat transfered in cell:\n')
		f.write('#'+'cell centre height [m]   transfered heat [kW] \n')
		for row in Q_trans:	
			f.write('	'.join(map(str,row)) + '\n')
		f.write('\n')

		f.write('#'+'Pressure data of both sides:\n')
		f.write('#'+'node height [m]   Primary Pressure [MPa] 	Secondary Pressure [MPa]\n')
		for row in Parray_data:	
			f.write('	'.join(map(str,row)) + '\n')

		f.write('#'+'velocity data of both sides:\n')
		f.write('#'+'node height [m]   Primary velocity [m/s] 	Secondary velocity [m/s]\n')
		for row in u_data:	
			f.write('	'.join(map(str,row)) + '\n')

		f.write('#'+'density data of both sides:\n')
		f.write('#'+'node height [m]   Primary density [kg/m^3] 	Secondary density [kg/m^3]\n')
		for row in rho_data:	
			f.write('	'.join(map(str,row)) + '\n')
	f.close()

# check if transient are needed
	if mode == 2:
		print ('steady state done')
		print (Tarray_pri[0],Tarray_sec[0])
		print (Tarray_pri[n-1],Tarray_sec[n-1])
		sys.exit()
	
#=========================restart from steady state=================================
if mode == 0:
	print ('\033[1m'+'**** restore data from steady state file ****'+'\033[0m')
#	restore data from file
	height_array = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = 0 ,skip_footer = (6*n-1))
	data_pri = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = 4 ,skip_footer = (5*n-1))
	data_sec = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = (n+1+5) ,skip_footer = (4*n-1))
	Q_array = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = 2*n+7, skip_footer = 3*n )
	P_array = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = 3*n+9, skip_footer = 2*n)
	u_array = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = 4*n+12, skip_footer = n)
	rho_array = np.genfromtxt(steadyfile,dtype = str, comments = '#',skip_header = 5*n+14 )
#	read node hight data
	height = [float(i) for i in data_pri[:,0]]	
#	reconstructe pressure, enthalpy, temperature, vapor fraction and transfered heated arrays
	height_boil = float(height_array[0])
	height_sat  = float(height_array[1])

	Parray_pri = [float(i) for i in P_array[:,1]]
	Parray_sec = [float(i) for i in P_array[:,2]]

	Tarray_pri = [float(i) for i in data_pri[:,1]]
	Tarray_sec = [float(i) for i in data_sec[:,1]]
	
	h_pri = [float(i) for i in data_pri[:,2]]
	h_sec = [float(i) for i in data_sec[:,2]]

	u_pri = [float(i) for i in u_array[:,1]]
	u_sec = [float(i) for i in u_array[:,2]]

	rho_pri = [float(i) for i in rho_array[:,1]]
	rho_sec = [float(i) for i in rho_array[:,2]]

	x_f = [float(i) for i in data_sec[:,3]]
	Q_transfer = [float(i) for i in Q_array[:,1]]

	Tarray_pri = np.array(Tarray_pri).T
	Tarray_sec = np.array(Tarray_sec).T

	Parray_pri = np.array(Parray_pri).T
	Parray_sec = np.array(Parray_sec).T

	h_pri = np.array(h_pri).T
	h_sec = np.array(h_sec).T

	u_pri = np.array(u_pri).T
	u_sec = np.array(u_sec).T

	rho_pri = np.array(rho_pri).T
	rho_sec = np.array(rho_sec).T

	x_f = np.array(x_f).T
	Q_transfer = np.array(Q_transfer).T


#=========================transient calculation=================================
print ('\033[1m'+'****start transient calculation ****'+'\033[0m')
# define the length of each region
l_l = height_boil
l_lg = height_sat - height_boil
l_g  = l_tube - height_sat

# record the steady state boiling and statured height
height_boil_steady = height_boil
height_sat_steady = height_sat

# number of timesteps
nstep = int((time_end - time_start)/dt)
# reset start time 
t_curr = time_start

# calculate the cell volume of each side
V_pri = A_pri*l_cell
V_sec = A_sec*l_cell

	
# initial Parray for last time step(state state)
Parray_pri_last = Parray_pri * np.ones(n)
Parray_sec_last = Parray_sec * np.ones(n)
# initial Tarray for last time step(state state)
Tarray_pri_last = Tarray_pri * np.ones(n)
Tarray_sec_last = Tarray_sec * np.ones(n)
# initial mass flow array for both sides for last timestep
Garray_pri_last = G_pri * np.ones(n)
Garray_sec_last = G_sec * np.ones(n)
# initial h for last time step(state state)
h_pri_last = h_pri * np.ones(n)
h_sec_last = h_sec * np.ones(n)
# initial rho for last time step(state state)
rho_pri_last = rho_pri * np.ones(n)
rho_sec_last = rho_sec * np.ones(n)
# initial u for last time step(state state)
u_pri_last = u_pri * np.ones(n)
u_sec_last = u_sec * np.ones(n)
# initial x for last time step(state state)
x_pri_last = np.zeros(n)
x_sec_last = x_f * np.ones(n)
# initial q array for transferd heat
q_last = Q_transfer * np.ones(n-1)

# current time step variable array
# initial Parray for current time step
Parray_pri_curr = Parray_pri*np.ones(n)
Parray_sec_curr = Parray_sec*np.ones(n)
# initial Tarray for current time step
Tarray_pri_curr = Tarray_pri*np.ones(n)
Tarray_sec_curr = Tarray_sec*np.ones(n)
#initial mass flow array current time step
Garray_pri_curr = G_pri*np.ones(n)
Garray_sec_curr = G_sec*np.ones(n)
# initial h for current time step
h_pri_curr = h_pri*np.ones(n)
h_sec_curr = h_sec*np.ones(n)
# initial rho for current time step
rho_pri_curr = rho_pri*np.ones(n)
rho_sec_curr = rho_sec*np.ones(n)
# initial u for current time step
u_pri_curr = u_pri*np.ones(n)
u_sec_curr = u_sec*np.ones(n)
# initial x for curr time step
x_pri_curr = np.zeros(n)
x_sec_curr = x_f*np.ones(n)
# initial q array for transferd heat
q_curr = Q_transfer*np.ones(n-1)

#print (rho_sec_last)
#sys.exit()

with open('Transfered_heat.txt','a') as f:
	f.write ('time '+ 'Transfered Heat (MW)'+'\n')

#time step iteration
for t in range(nstep):
	t_curr = t_curr + dt
	if t_curr > time_end:
		print ('***************')
		print ('Simulation Done')
		print ('***************')
		break
	print ('************************************')
	print ('current time', t_curr)
	print ('************************************')

#	input information of primary side
#	T_inlet_pri = T_inlet_pri + T_ck
#	h_inlet_pri = IAPWS97(P = Parray_pri_curr[0],T = T_inlet_pri).h

#	update boundary condition, primary side
	Tarray_pri_curr[n-1] = T_inlet_pri
	Parray_pri_curr[n-1] = P_pri
	Garray_pri_curr[n-1] = Garray_pri_last[n-1]
#	calcualte boundary state variable
	rho_pri_curr[n-1] = IAPWS97(P = Parray_pri_curr[n-1], T = Tarray_pri_curr[n-1]).rho
	u_pri_curr[n-1] = Garray_pri_curr[n-1]/(rho_pri_curr[n-1]*A_pri)
	h_pri_curr[n-1] = IAPWS97(P = Parray_pri_curr[n-1], T = Tarray_pri_curr[n-1]).h 
	 
#	update boundary condition, secondary side
	Tarray_sec_curr[0] = Tarray_sec_last[0] 
	Parray_sec_curr[0] = max(Parray_sec_last[0] - 0.05/(1/dt) , Parray_sec[0] - 3.6)
	Garray_sec_curr[0] = Garray_sec_last[0] 
#	if abs(Parray_sec_last[0] - Parray_sec[0]) >= 0.2:
#		print ('Abnormal')
#		try:
#			t_trigger
#		except NameError:
#			t_trigger = t_curr + t_iso + t_trip 
#			print ('Isolation start')
#		if t_curr - t_trigger >= 0:
#			print ('Isolating')
#			Garray_sec_curr[0] = max(G_sec*math.exp(-0.5*(t_curr - t_trigger)) , 0.01*G_sec)
#		else:
#			Garray_sec_curr[0] = Garray_sec_last[0] 
#	else:
#		Garray_sec_curr[0] = Garray_sec_last[0] 
#	calcualte boundary state variable
	rho_sec_curr[0] ='%.10f'%IAPWS97(P = Parray_sec_curr[0], T = Tarray_sec_curr[0]).rho
	u_sec_curr[0] = '%.10f'%(Garray_sec_curr[0]/(rho_sec_curr[0]*A_sec))
	h_sec_curr[0] = '%.10f'%IAPWS97(P = Parray_sec_curr[0], T = Tarray_sec_curr[0]).h

#	reset the phase condition
	h_boil = float(format(IAPWS97(P = Parray_sec_curr[0], x = 0).h, '.10f'))
	h_sat = float(format(IAPWS97(P = Parray_sec_curr[0], x = 1).h, '.10f'))


#	initial heat transfer
	q_curr = q_last * np.ones(n-1)
	q_temp = q_last * np.ones(n-1)

	for j in range(n_loop):	

	#	data primary side
		for i in range(n-1):
			u_pri_old = u_pri_curr[n-i-1-1] 
			rho_pri_old = rho_pri_curr[n-i-1-1] 
			h_pri_old = h_pri_curr[n-i-1-1]
			for k in range(n_loop):
			#	cell average data
				T_pri_m_curr = 1/2 * (Tarray_pri_curr[n-i-1] + Tarray_pri_curr[n-i-1-1])
				P_pri_m_curr = 1/2 * (Parray_pri_curr[n-i-1] + Parray_pri_curr[n-i-1-1])
				G_pri_m_curr = 1/2 * (Garray_pri_curr[n-i-1] + Garray_pri_curr[n-i-1-1])
				x_pri_m_curr = 1/2 * (x_pri_curr[n-i-1] + x_pri_curr[n-i-1-1])

			#	total pressure drop by corrolations
				Parray_pri_curr[n-i-1-1] = Parray_pri_curr[n-i-1] - pd.pd_single(Garray_pri_curr[n-i-1-1],Garray_pri_curr[n-i-1],T_pri_m_curr,P_pri_m_curr,rho_pri_curr[n-i-1-1],rho_pri_curr[n-i-1],x_pri_m_curr,dh_pri,A_pri,l_cell,-1)
		
			#	slove energy equation
				h_pri_curr[n-i-1-1] = eq.h_curr(h_pri_curr[n-i-1],h_pri_last[n-i-1-1],h_pri_last[n-i-1],u_pri_curr[n-i-1-1],u_pri_curr[n-i-1],rho_pri_curr[n-i-1-1],rho_pri_curr[n-i-1],l_cell,dt,A_pri,-q_curr[n-i-1-1])

			#	lookup data from IAPWS97
				Tarray_pri_curr[n-i-1-1] = float('%.10f'%IAPWS97( P = Parray_pri_curr[n-i-1-1], h = h_pri_curr[n-i-1-1]).T)
				rho_pri_curr[n-i-1-1] = float('%.10f'%IAPWS97(P = Parray_pri_curr[n-i-1-1], h = h_pri_curr[n-i-1-1]).rho)
		
			#	solve momentum equation 
				u_pri_curr[n-i-1-1] = eq.u_correct(u_pri_curr[n-i-1],rho_pri_curr[n-i-1-1],rho_pri_curr[n-i-1],rho_pri_last[n-i-1-1],rho_pri_curr[n-i-1],dt,l_cell)
			#	u_pri_curr[n-i-1-1] = float('%.10f'%(G_pri/(rho_pri_curr[n-i-1]*A_pri)))
		
			#	calculate new mass flow rate
				Garray_pri_curr[n-i-1-1] = float('%.10f'%(rho_pri_curr[n-i-1-1] * u_pri_curr[n-i-1-1] * A_pri))
	
				if (abs(rho_pri_old - rho_pri_curr[n-i-1-1]) < eps_trans*10 and abs(u_pri_old - u_pri_curr[n-i-1-1]) < eps_trans
					and abs(h_pri_old - h_pri_curr[n-i-1-1] < eps_trans*10)):
					break
				elif k >= (n_loop-1):
					print ('primary unconverge')
					sys.exit()
				rho_pri_old = rho_pri_curr[n-i-1-1]
				u_pri_old = u_pri_curr[n-i-1-1]
				h_pri_old = h_pri_curr[n-i-1-1]
#		print (Tarray_pri_curr - Tarray_pri_last)
#		print (Parray_pri_curr - Parray_pri_last)
#		sys.exit()	

	#	data secondary side
		for i in range(n-1):
			u_sec_old = u_sec_curr[i+1] 
			rho_sec_old = rho_sec_curr[i+1]
			h_sec_old = h_sec_curr[i+1]
			q_old = q_curr[i]
		
			for k in range(n_loop):
			#	cell average temperature for primary side
				T_pri_m_curr = 1/2 * (Tarray_pri_curr[i] + Tarray_pri_curr[i+1])
				P_pri_m_curr = 1/2 * (Parray_pri_curr[i] + Parray_pri_curr[i+1])
				G_pri_m_curr = 1/2 * (Garray_pri_curr[i] + Garray_pri_curr[i+1])
				x_pri_m_curr = 1/2 * (x_pri_curr[i] + x_pri_curr[i+1])
			#	cell average data for secondary side
				T_sec_m_curr = 1/2 * (Tarray_sec_curr[i] + Tarray_sec_curr[i+1])
				P_sec_m_curr = 1/2 * (Parray_sec_curr[i] + Parray_sec_curr[i+1])
				G_sec_m_curr = 1/2 * (Garray_sec_curr[i] + Garray_sec_curr[i+1])
				x_sec_m_curr = 1/2 * (x_sec_curr[i] + x_sec_curr[i+1])

				if x_sec_m_curr == 0 or x_sec_m_curr >= x_b:
				#	single phase pressure drop
					Parray_sec_curr[i+1] = Parray_sec_curr[i] - pd.pd_single(Garray_sec_curr[i+1],Garray_sec_curr[i],T_sec_m_curr,P_sec_m_curr,rho_sec_curr[i+1],rho_sec_curr[i],x_sec_m_curr,dh_sec,A_sec,l_cell,1)
				elif x_sec_m_curr > 0 and x_sec_m_curr < x_b:
				#	two phase pressure drop
					Parray_sec_curr[i+1] = Parray_sec_curr[i] - pd.pd_twHEM(Garray_sec_curr[i+1],Garray_sec_curr[i],T_sec_m_curr,P_sec_m_curr,rho_sec_curr[i+1],rho_sec_curr[i],x_sec_curr[i+1],x_sec_curr[i],dh_sec,A_sec,l_cell,1)
	
			#	calculate heat source or sink for the cell
			#	log mean temperature
				dT_up = Tarray_pri_curr[i+1] - Tarray_sec_curr[i+1]
				dT_low = Tarray_pri_curr[i] - Tarray_sec_curr[i]
				dT_m = (dT_up - dT_low)/math.log(dT_up/dT_low)
			#	calcualte secondary wall temperature
				Twall = hc.Tw(G_pri_m_curr,A_pri,r_pri,r_sec,l_cell,dh_pri,n_tube,T_pri_m_curr,P_pri_m_curr,q_curr[i])
			#	calculate overall heat transfer coefficient
				UA_pri = hc.UA(G_pri_m_curr,A_pri,r_pri,r_sec,l_cell,dh_pri,T_pri_m_curr,P_pri_m_curr,x_pri_curr[i+1],x_pri_curr[i], n_tube, q_curr[i],Twall[0])
				UA_sec = hc.UA(G_sec_m_curr,A_sec,r_sec,r_pri,l_cell,dh_sec,T_sec_m_curr,P_sec_m_curr,x_sec_curr[i+1],x_sec_curr[i], n_tube, q_curr[i],Twall[1])
				UA = 1/(1/UA_pri + 1/UA_sec)
			#	waring for nonsence value	
				if UA <= 0:
					print ('************* WARNING **************')
					print ('**negative heat transfer coeffient**')
					print ('************************************')
			#	calculate transfered heat#		
				q = float('%.10f'%(UA*dT_m/1000 * n_tube)) # in unit kW
				q_curr[i] = q
#				print (q_curr[i] - q_last[i])
#				sys.exit()

			#	slove momentum equation, mass continuty and energy equation
				h_sec_curr[i+1] = eq.h_curr(h_sec_curr[i],h_sec_last[i+1],h_sec_last[i],u_sec_curr[i+1],u_sec_curr[i],rho_sec_curr[i+1],rho_sec_curr[i],l_cell,dt,A_sec,q_curr[i])

			#	lookup data from IAPWS97	
				Tarray_sec_curr[i+1] = float('%.10f'%IAPWS97( P = Parray_sec_curr[i+1], h = h_sec_curr[i+1]).T)
				rho_sec_curr[i+1] = float('%.10f'%IAPWS97( P = Parray_sec_curr[i+1], h = h_sec_curr[i+1]).rho)
				x_sec_curr[i+1] = float('%.10f'%(IAPWS97( P = Parray_sec_curr[i+1], h = h_sec_curr[i+1]).x))
	
			#	solve mass continuity equation for velocity	
				u_sec_curr[i+1] = eq.u_correct(u_sec_curr[i],rho_sec_curr[i+1],rho_sec_curr[i],rho_sec_last[i+1],rho_sec_curr[i],dt,l_cell)
			#	u_sec_curr[i+1] = float('%.10f'%(G_sec/(rho_sec_curr[i+1]*A_sec)))
			#	calculate mass flow rate
				Garray_sec_curr[i+1] = float('%.10f'%(rho_sec_curr[i+1] * u_sec_curr[i+1] * A_sec))
	
			#	recalculate the boiling enthalpy
				h_boil = float('%.10f'%(IAPWS97(P = Parray_sec[i+1], x = 0).h))
			#	recalculate the saturated enthalpy
				h_sat = float('%.10f'%(IAPWS97(P = Parray_sec[i+1], x = 1).h))
	
				if (abs(rho_sec_old - rho_sec_curr[i+1]) < eps_trans*10 and  abs(u_sec_old - u_sec_curr[i+1]) < eps_trans
						and abs(h_sec_old - h_sec_curr[i+1]) < eps*10 and abs(q_old - q) < eps*100):
					q_curr[i] = q
					break
				elif k >= (n_loop-1):
					print ('seccondary unconverge')
					sys.exit()
				rho_sec_old = rho_sec_curr[i+1]
				u_sec_old = u_sec_curr[i+1]
				h_sec_old = h_sec_curr[i+1]
				q_old = q
#		print (Tarray_sec_curr - Tarray_sec_last)
#		print (Parray_sec_curr - Parray_sec_last)
#		print (Garray_sec_curr - Garray_sec_last)
	#	print (q_curr)
		Q_trans = sum(q_curr)
		Q_pri = (h_pri_curr[n-1] - h_pri_curr[0])*(G_pri/1000)
		Q_sec = (h_sec_curr[n-1] - h_sec_curr[0])*(G_sec/1000)
		
		print (Q_trans/1000,Q_pri,Q_sec)
		print (Tarray_pri_curr[n-1],Tarray_sec_curr[n-1])
		print (Tarray_pri_curr[0],Tarray_sec_curr[0])

		con_check = abs(q_curr - q_temp)
		if all (item < eps_trans for item in con_check):
			print ('number of iterations', j)
	#		print (Tarray_pri_curr - Tarray_pri)
		#	print (Parray_pri_curr - Parray_pri_last)
		#	print (Garray_pri_curr - Garray_pri_last)
		#	print (Tarray_sec_curr - Tarray_sec)
		#	print (Parray_sec_curr - Parray_sec_last)
		#	print (Garray_sec_curr - Garray_sec_last)
		#	print (q_curr - q_temp)
			break	
		elif j >= (n-1):
			print ('Heat transfer unconverge')
			break
		q_temp = q_curr * np.ones(n-1)

#	for i in range(n-1):
#		print('time:','%.4f'%t_curr,'node number:' ,i,'height:','%.10f'%((i)*l_cell))
#		print ('Primary absolute pressure',Parray_pri[i],Parray_pri_curr[i]) 
#		print ('Secondary absolute pressure',Parray_sec[i],Parray_sec_curr[i]) 
#		print('primary temperature',Tarray_pri[i],Tarray_pri_curr[i])
#		print('secdonary temperature',Tarray_sec[i],Tarray_sec_curr[i])
#		print('Primary mass flow rate',G_pri,Garray_pri_curr[i])
#		print('secdonary mass flow rate',G_sec,Garray_sec_curr[i])
#		print('secdonary vapor fraction',x_f[i],x_sec_curr[i])
#		print('primary velocity',u_pri[i],u_pri_curr[i])
#		print('secdonary velocity',u_sec[i],u_sec_curr[i])
#		print('primary density',rho_pri[i],rho_pri_curr[i])
#		print('secdonary density',rho_sec[i],rho_sec_curr[i])
#		print ('\n')
	#	determine boiling
#		if x_sec_curr[i+1] > 0 and x_sec_curr[i] == 0:
	#		node of boiling and height of boilling
#			height_boil = (i+1)*l_cell - l_cell/2
#			print ('node boiling', i+1, 'at height', '%.5f'%height_boil, '[m]')
#			print ('steady boiling at height','%.5f'%height_boil_steady, '[m]')
#	
#	#	determine saturate steam
#		elif x_sec_curr[i] < x_b and x_sec_curr[i+1] >= x_b:
#			height_sat = (i+1)*l_cell - l_cell/2
#			print ('node saturated', i+1, 'at height','%.5f'%height_sat, '[m]')
#			print ('steady saturated at height','%.5f'%height_sat_steady, '[m]')
#		print ('\n')

	for i in range(n-1):
	#	determine boiling
		if x_sec_curr[i+1] > 0 and x_sec_curr[i] == 0:
	#		node of boiling and height of boilling
			height_boil = (i+1)*l_cell - l_cell/2
			print ('node boiling', i+1, 'at height', '%.10f'%height_boil, '[m]')
	
	#	determine saturate steam
		elif x_sec_curr[i] < x_b and x_sec_curr[i+1] >= x_b:
			height_sat = (i+1)*l_cell - l_cell/2
			print ('node saturated', i+1, 'at height','%.10f'%height_sat, '[m]')

	print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')	
	print('total transfered heat', Q_trans/1000, "MW")
	print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')	
	print('\n')

	# calculate temperature difference
	dT_pri = Tarray_pri_curr - Tarray_pri
	dT_sec = Tarray_sec_curr - Tarray_sec
	# calculate Mass flow difference
	dG_pri = Garray_pri_curr - G_pri
	dG_sec = Garray_sec_curr - G_sec
	# calculate Pressure difference
	dP_pri = Parray_pri_curr - Parray_pri
	dP_sec = Parray_sec_curr - Parray_sec
	# calculate enthalpy difference
	deltah_pri = h_pri_curr - h_pri
	deltah_sec = h_sec_curr - h_sec
	# calculate density difference
	drho_pri = rho_pri_curr - rho_pri
	drho_sec = rho_sec_curr - rho_sec
	# calculate density difference
	du_pri = u_pri_curr - u_pri
	du_sec = u_sec_curr - u_sec

	with open(transfile, 'a') as f:
		f.write ('*****************'+'\n')
		f.write ('time '+ str(t_curr)+'\n')
		f.write ('*****************'+'\n')
		f.write('Height   Primary Temperature   Secondary Temperature'+'\n')
		data = np.column_stack((height,Tarray_pri_curr,Parray_pri_curr,Tarray_sec_curr,Parray_sec_curr))
		for row in data:	
			f.write('	'.join(map(str,row)) + '\n')
	with open('Massflow_compare.txt','a') as f:
		f.write ('*****************'+'\n')
		f.write ('time '+ str(t_curr)+'\n')
		f.write ('*****************'+'\n')
		data = np.column_stack((height,Garray_pri_curr,Garray_sec_curr,u_pri_curr,u_sec_curr))
		for row in data:	
			f.write('	'.join(map(str,row)) + '\n')
	with open('enthalpy_compare.txt','a') as f:
		f.write ('*****************'+'\n')
		f.write ('time '+ str(t_curr)+'\n')
		f.write ('*****************'+'\n')
		data = np.column_stack((height,h_pri_curr,h_pri,h_sec_curr,h_sec))
		for row in data:	
			f.write('	'.join(map(str,row)) + '\n')
	with open('velocity_compare.txt','a') as f:
		f.write ('*****************'+'\n')
		f.write ('time '+ str(t_curr)+'\n')
		f.write ('*****************'+'\n')
		data = np.column_stack((height,u_pri_curr,u_pri,u_sec_curr,u_sec))
		for row in data:	
			f.write('	'.join(map(str,row)) + '\n')
	with open('Transfered_heat.txt','a') as f:
		Q_pri = (h_pri_curr[n-1] - h_pri_curr[0])*(G_pri/1000)
		Q_sec = (h_sec_curr[n-1] - h_sec_curr[0])*(G_sec/1000)
		f.write('%.5f'%t_curr+'  '+'%.5f'%Q_pri+'  '+'%.5f'%Q_sec+'  '+'%.5f'%(Q_trans/1000)+'\n')
	with open('InoutT.txt','a') as f:
		f.write('%.5f'%t_curr+'  '+'%.3f'%Tarray_pri_curr[0]+'  '+'%.3f'%Tarray_pri_curr[n-1]+'  '+'%.3f'%Tarray_sec_curr[0]+'  '+'%.3f'%Tarray_sec_curr[n-1]+'\n')
	with open('InoutP.txt','a') as f:
		f.write('%.5f'%t_curr+'  '+'%.5f'%Parray_pri_curr[0]+'  '+'%.5f'%Parray_pri_curr[n-1]+'  '+'%.5f'%Parray_sec_curr[0]+'  '+'%.5f'%Parray_sec_curr[n-1]+'\n')
	with open('Phasechange_Loc.txt','a') as f:
		f.write('%.5f'%t_curr+'  '+'%.5f'%height_boil+'  '+'%.5f'%height_sat+'\n')
	with open('udiff.txt','a') as f:
		f.write('%.5f'%t_curr+'  '+'%.5f'%du_pri[0]+'  '+'%.5f'%du_pri[n-1]+'  '+'%.5f'%du_sec[0]+'  '+'%.3f'%du_sec[n-1]+'\n')

	# mass flow array for both sides for last timestep
	Garray_pri_last = Garray_pri_curr * np.ones(n)
	Garray_sec_last = Garray_sec_curr * np.ones(n)
	
	# Parray for last time step
	Parray_pri_last = Parray_pri_curr * np.ones(n)
	Parray_sec_last = Parray_sec_curr * np.ones(n)
	
	# Tarray for last time step
	Tarray_pri_last = Tarray_pri_curr * np.ones(n)
	Tarray_sec_last = Tarray_sec_curr * np.ones(n)
	
	# h for last time step
	h_pri_last = h_pri_curr * np.ones(n)
	h_sec_last = h_sec_curr * np.ones(n)
	
	# rho for last time step
	rho_pri_last = rho_pri_curr * np.ones(n)
	rho_sec_last = rho_sec_curr * np.ones(n)
	
	# u for last time step
	u_pri_last = u_pri_curr * np.ones(n)
	u_sec_last = u_sec_curr * np.ones(n)
	
	# x for last time step
	x_sec_last = x_sec_curr * np.ones(n)
	x_sec_last = x_sec_curr * np.ones(n)
	
	# transfered heat for last time step
	q_last = q_curr*np.ones(n-1)
