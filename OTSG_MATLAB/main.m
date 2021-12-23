H_SG = 0:0.1:6.1;
PD_Ratio = 1.3;
P_Sec = 40;
T_Sec_feed = 125.85;
T_Pri_Set(1,1:62) = 0;
T_Sec_Set(1,1:62) = 0;
for mmmmm = 0:100
    T_Sec_feed = T_Sec_feed + 1;
    TCalculation;
    T_Pri_Set(mmmmm+1,:) = T_Pri;
    T_Sec_Set(mmmmm+1,:) = T_Sec;
    clear T_Pri T_Sec
end