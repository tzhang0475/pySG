%H_SG = 0:0.1:6.1;
%PD_Ratio = 1.3;
Power = 2e8;
P_Pri = 145;
%P_Sec = 40;
%T_Sec_feed = 120;
T_Pri(1) = 271;
T_Sec(1) = T_Sec_feed;
q = 0;
Saturation = 0;
m_dot_Pri = 1000;
m_dot_Sec = Power / 1000 / (XSteam('h_pT',P_Sec,XSteam('Tsat_p',P_Sec)+30)-XSteam('h_pT',P_Sec,T_Sec_feed));
%Calculation of Liquid Phase:
for mm = 1:61
    T_Pri(mm+1) = T_Pri(mm) + 0.6;
    T_Sec(mm+1) = T_Sec(mm) + 3;
    if T_Sec(mm+1) >= XSteam('Tsat_P',P_Sec)
        Saturation = 1;
        break
    end
    for nn = 1:1000
        T_Pri_m = 0.5 * ( T_Pri(mm) + T_Pri(mm+1));
        T_Sec_m = 0.5 * ( T_Sec(mm) + T_Sec(mm+1));
        [UA,m_dot_Pri,m_dot_Sec] = f_LHTC(PD_Ratio,T_Pri_m,P_Sec,T_Sec_m,T_Sec_feed);
        delta_T1 = T_Pri(mm+1) - T_Sec(mm+1);
        delta_T2 = T_Pri(mm) - T_Sec(mm);
        delta_Tlm = (delta_T1-delta_T2)/log(delta_T1/delta_T2);
        q = UA * delta_Tlm/1000*4237;
        T_Pri_new = XSteam('T_ph',P_Pri,q/m_dot_Pri+XSteam('h_pT',P_Pri,T_Pri(mm)));
        T_Sec_new = XSteam('T_ph',P_Sec,q/m_dot_Sec+XSteam('h_pT',P_Sec,T_Sec(mm)));
        if T_Sec_new >= XSteam('Tsat_P',P_Sec)
            Saturation = 1;
            break
        end
        if abs(T_Pri_new -T_Pri(mm+1))<=1e-8 && abs(T_Sec_new - T_Sec(mm+1))<=1e-8
            T_Pri(mm+1) = T_Pri_new;
            T_Sec(mm+1) = T_Sec_new;
            break
        else
            T_Pri(mm+1) = T_Pri_new;
            T_Sec(mm+1) = T_Sec_new;
        end
    end
    if(Saturation==1)
        break
    end
end
%Assume saturation water node:
T_Sec(mm+1) = XSteam('Tsat_p',P_Sec);
T_Pri(mm+1) = XSteam('T_ph',P_Pri,XSteam('h_pT',P_Pri,T_Pri(mm)) ...
    +(XSteam('hL_p',P_Sec)-XSteam('h_pT',P_Sec,T_Sec(mm)))*m_dot_Sec/m_dot_Pri);
mm = mm + 1;
h_Sec(mm) = XSteam('hL_p',P_Sec);
X_Sec(mm) = 0;
%Calculation of Boiling:
for mmm = mm:61
    q = 2000; %kW
    for nn = 1:1000
        T_Pri(mmm+1) = XSteam('T_ph',P_Pri,XSteam('h_pT',P_Pri,T_Pri(mmm))+q/m_dot_Pri);
        T_Pri_m = 0.5 * ( T_Pri(mmm) + T_Pri(mmm+1));
        T_Sec(mmm+1) = T_Sec(mmm);
        T_Sec_m = 0.5 * (T_Sec(mmm) + T_Sec(mmm+1));
        h_Sec(mmm+1) = h_Sec(mmm) + q/m_dot_Sec;
        X_Sec(mmm+1) = XSteam('x_ph',P_Sec,h_Sec(mmm+1));
        if X_Sec(mmm+1) == 1
            Saturation = 2;
            break
        end
        X_Sec_m = 0.5 * (X_Sec(mmm) + X_Sec(mmm+1));
        [UA,m_dot_Pri,m_dot_Sec] = f_BHTC(PD_Ratio,T_Pri_m,P_Sec,T_Sec_m,X_Sec_m,T_Sec_feed,q);
        delta_T1 = T_Pri(mmm+1) - T_Sec(mmm+1);
        delta_T2 = T_Pri(mmm) - T_Sec(mmm);
        delta_Tlm = (delta_T1-delta_T2)/log(delta_T1/delta_T2);
        q_new = UA * delta_Tlm/1000*4237;
        if abs(q - q_new) <= 1e-8
            q = q_new;
            break
        else
            q = q_new;
        end
    end
    if Saturation == 2
        break
    end
end
%Assume saturation steam node:
T_Sec(mmm+1) = XSteam('Tsat_p',P_Sec);
T_Pri(mmm+1) = XSteam('T_ph',P_Pri,XSteam('h_pT',P_Pri,T_Pri(mmm)) ...
    +(XSteam('hV_p',P_Sec)-XSteam('h_px',P_Sec,X_Sec(mmm)))*m_dot_Sec/m_dot_Pri);
mmm = mmm + 1;
h_Sec(mmm) = XSteam('hV_p',P_Sec);
X_Sec(mmm) = 1;
%Calculation of Overheating:
for mmmm = mmm:61
    T_Pri(mmmm+1) = T_Pri(mmmm) + 0.3;
    T_Sec(mmmm+1) = T_Sec(mmmm) + 2.0;
    for nn = 1:1000
        T_Pri_m = 0.5 * ( T_Pri(mmmm) + T_Pri(mmmm+1));
        T_Sec_m = 0.5 * ( T_Sec(mmmm) + T_Sec(mmmm+1));
        [UA,m_dot_Pri,m_dot_Sec] = f_VHTC(PD_Ratio,T_Pri_m,P_Sec,T_Sec_m,T_Sec_feed);
        delta_T1 = T_Pri(mmmm+1) - T_Sec(mmmm+1);
        delta_T2 = T_Pri(mmmm) - T_Sec(mmmm);
        delta_Tlm = (delta_T1-delta_T2)/log(delta_T1/delta_T2);
        q = UA * delta_Tlm/1000*4237;
        T_Pri_new = XSteam('T_ph',P_Pri,q/m_dot_Pri+XSteam('h_pT',P_Pri,T_Pri(mmmm)));
        T_Sec_new = XSteam('T_ph',P_Sec,q/m_dot_Sec+h_Sec(mmmm));
        if abs(T_Pri_new -T_Pri(mmmm+1))<=1e-8 && abs(T_Sec_new - T_Sec(mmmm+1))<=1e-8
            T_Pri(mmmm+1) = T_Pri_new;
            T_Sec(mmmm+1) = T_Sec_new;
            h_Sec(mmmm+1) = XSteam('h_pT',P_Sec,T_Sec(mmmm+1));
            break
        else
            T_Pri(mmmm+1) = T_Pri_new;
            T_Sec(mmmm+1) = T_Sec_new;
            h_Sec(mmmm+1) = XSteam('h_pT',P_Sec,T_Sec(mmmm+1));
        end
    end
end