# solver for continuty, momentum and energy equation

# libs for calculation 
import numpy as np
import pdrop as pd
import math


# steady state momentum equation
def u_steady(u_i,rho_ip1,rho_i,dp):
	rho_m = 1/2*(rho_ip1 + rho_i)
	result = '%.10f'%(math.sqrt(u_i**2 - 2/rho_m * dp))

	return result

# steady state energy equation
def rho_steady(rho_i,u_ip1,u_i):
	u_m = 1/2 * (u_ip1+u_i)
	du = u_ip1 - u_i

	a = du/2 + u_m
	c = (u_m - du/2)*rho_i

	result = c/a
	
	return '%.10f'%result

# solution of rho for continuty equation
def rho_curr(u_n_ip1,u_n_i,rho_n_i,rho_nm1_ip1,rho_nm1_i,dz,dt):
	u_n_m = 1/2 * (u_n_ip1 + u_n_i)
	rho_nm1_m = 1/2 * (rho_nm1_ip1+rho_nm1_i)
	du_n_z = u_n_ip1 - u_n_i
	b = (dt*u_n_m-dz/2-dt*du_n_z/2)*rho_n_i + dz*rho_nm1_m 
	a = dz/2 + dt*du_n_z/2 + dt*u_n_m
	result = '%.10f'%(b/a)
	
	return result

# solution of u for mass continuty equation
def u_correct(u_n_i,rho_n_ip1,rho_n_i,rho_nm1_ip1,rho_nm1_i,dt,dz):
	rho_n_m = 1/2 * (rho_n_ip1 + rho_n_i)
	rho_nm1_m = 1/2 * (rho_nm1_ip1 + rho_nm1_i)
	drho_n_z = rho_n_ip1 - rho_n_i
	drho_t = rho_n_m - rho_nm1_m

	b1 =  (rho_n_m - 1/2*drho_n_z)*u_n_i
	b2 = -dz*drho_t/dt
	a = rho_n_m + 1/2*drho_n_z

	result = '%.10f'%((b1+b2)/a)

	return result

# solution of velocity for momentum equation
def u_curr(u_n_i,u_nm1_ip1,u_nm1_i,rho_n_ip1,rho_n_i,dz,dt,dp):
#	average filed
	u_nm1_m = 1/2 * (u_nm1_ip1 + u_nm1_i)
	rho_n_m = 1/2 * (rho_n_ip1 + rho_n_i)
#	parameter for equation
	a = dt
	b = dz
	c1 = -dt*u_n_i**2 + dz*u_n_i - 2*dz*u_nm1_m
	c2 = 2*dt*dp/rho_n_m
	c = c1+c2

#	general solution for quadratic equation
	slv_1 = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)
	slv_2 = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)


	result = '%.10f'%max(slv_1,slv_2)

	return result
# solution of enalthpy for single phase energy equation
def h_single(h_n_i,h_nm1_ip1,h_nm1_i,u_n_ip1,u_n_i,rho_n_ip1,rho_n_i,dz,dt,q):
	a = dt*q/(1/2*(rho_n_ip1 + rho_n_i))
	b = dz*(1/2)*(h_nm1_ip1+h_nm1_i) + dt*(1/2)*(u_n_ip1+u_n_i)*h_n_i - dz*(1/2)*h_n_i
	c = dz/2 + dt*(1/2)*(u_n_ip1+u_n_i) 

	result = '%.10f'%((a + b)/c)
	
	return result

# solution of enalthpy for two two phase energy equation
def h_tw(h_n_i,h_nm1_ip1,h_nm1_i,u_n_ip1,u_n_i,rho_n_ip1,rho_n_i,dz,dt,q,A,dp):
#	a = (dt*q +dp*A*dz*(dt+dz))/(1/2*(rho_n_ip1 + rho_n_i))  
	a = dt*q/(1/2*(rho_n_ip1 + rho_n_i))
	b = dz*(1/2)*(h_nm1_ip1+h_nm1_i) + dt*(1/2)*(u_n_ip1+u_n_i)*h_n_i - dz*(1/2)*h_n_i
	c = dz/2 + dt*(1/2)*(u_n_ip1+u_n_i) 

	result = '%.10f'%((a + b)/c)

	return result

# solution of enalthpy for energy equation
def h_curr(h_n_i,h_nm1_ip1,h_nm1_i,u_n_ip1,u_n_i,rho_n_ip1,rho_n_i,dz,dt,A_xsc,q):
#	volume heat source
	q_v = q/(A_xsc*dz)
#	average field
	u_n_m = 1/2 * (u_n_ip1 + u_n_i)
	rho_n_m = 1/2 * (rho_n_ip1 + rho_n_i)
	h_nm1_m = 1/2 * (h_nm1_ip1 + h_nm1_i)
#	calculate parameters
	a = dz/2 + dt*u_n_m
	b = dz*h_nm1_m + (dt*u_n_m - dz/2)*h_n_i 
	c = dt*dz*q_v/rho_n_m

	result = '%.10f'%((b+c)/a)
	
	return result

