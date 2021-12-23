function Nu = f_Nu_Gnielinski(Re,Pr)
f = f_f_Colebrook(Re);
Nu = (f/8)*(Re-1000)*Pr/(1+12.7*(f/8)^(0.5)*(Pr^(2/3)-1));