function [Gc_Sys,H_Sys,Kp,Ti,Td]=opt_pid(key,typ,vars)
%
% Optimum PID Controller Design &  ultimate frequency and gain where
% key = 2, 3, 4 for PI, normal PID, and PID controllers with D in the feedback path(only for set point) and
% typ = 1, 2 for set-point and disturbance rejection, respectively. The variable 
% vars = [k, L, T, N, C] for Optimum PID where C is the criterion type with C = 1, 2, 3
%        for ISE, ISTE, and IST2E criteria, respectively
%      & [Kc,Tc,kappa] for ultimate frequency and gain based on ISTE criteria.
%
% Disturbance Rejection part has been added to this function by: 
% Sina Masoumi, BSc student, Dept of Mech Eng., Sharif U of T., Dec. 2020, Dey 1399.
%
Ti=0; Td=0; N=vars(4); 
if length(vars)==5, k=vars(1); L=vars(2); T=vars(3); iC=vars(5); tt=0;
else, Kc=vars(1); Tc=vars(2); kappa=(3); tt=1; end
if tt==0
    if typ == 1
           switch key
           case 2
              PIDtab=[0.980, 0.712, 0.569, 1.072, 0.786, 0.628;
                     -0.892,-0.921,-0.951,-0.560,-0.559,-0.583;
                      0.690, 0.968, 1.023, 0.648, 0.883, 1.007;
                     -0.155,-0.247,-0.179,-0.114,-0.158,-0.167];
           case 3
              PIDtab=[1.048, 1.042, 0.968, 1.154, 1.142, 1.061; 
                     -0.897,-0.897,-0.904,-0.567,-0.579,-0.583;
                      1.195, 0.987, 0.977, 1.047, 0.919, 0.892;
                     -0.368,-0.238,-0.253,-0.220,-0.172,-0.165;
                      0.489, 0.385, 0.316, 0.490, 0.384, 0.315;
                      0.888, 0.906, 0.892, 0.708, 0.839, 0.832];
           case 4
              PIDtab=[1.260, 1.053, 0.942, 1.295, 1.120, 1.001;
                     -0.887,-0.930,-0.933,-0.619,-0.625,-0.624;
                      0.701, 0.736, 0.770, 0.661, 0.720, 0.754;
                     -0.147,-0.126,-0.130,-0.110,-0.114,-0.116;
                      0.375, 0.349, 0.308, 0.378, 0.350, 0.308;
                      0.886, 0.907, 0.897, 0.756, 0.811, 0.813];
           end
           ii=0; 
           if (L/T>1) 
               ii=3; 
           end 
           tt=L/T; 
           a1=PIDtab(1,ii+iC); b1=PIDtab(2,ii+iC); a2=PIDtab(3,ii+iC); b2=PIDtab(4,ii+iC); 
           Kp=a1/k*tt^b1; Ti=T/(a2+b2*tt); 
           if key==3|| key==4
               a3=PIDtab(5,ii+iC); b3=PIDtab(6,ii+iC); Td=a3*T*tt^b3; 
           end
    elseif typ == 2 
        switch key
            case 2
                PIDtab=[1.279 1.015 1.021 1.346 1.065 1.076;
                -0.945 -0.957 -0.953 -0.675 -0.673 -0.648;
                0.535 0.667 0.629 0.552 0.687 0.650;
                0.586 0.552 0.546 0.438 0.427 0.442];
           case 3
                 PIDtab=[1.473 1.468 1.531 1.524 1.515 1.592;
                -0.970 -0.970 -0.960 -0.735 -0.730 -0.705;
                1.115 0.942 0.971 1.130 0.957 0.957;
                0.753 0.725 0.746 0.641 0.598 0.597;
                0.550 0.443 0.413 0.552 0.444 0.414;
                0.948 0.939 0.933 0.851 0.847 0.850];
         end
           ii=0; 
           if (L/T>1) 
               ii=3; 
           end 
           tt=L/T; 
           a1=PIDtab(1,ii+iC); b1=PIDtab(2,ii+iC); a2=PIDtab(3,ii+iC); b2=PIDtab(4,ii+iC); 
           Kp=a1/T*tt^b1; Ti=T/a2*tt^b2; 
           if key==3
               a3=PIDtab(5,ii+iC); b3=PIDtab(6,ii+iC); Td=a3*T*tt^b3; 
           end
    end
elseif tt == 1
    if typ == 1
        switch key
            case 2,Kp=(4.264-0.184*kappa)/(12.119-0.432*kappa)*Kc...
                    ; Ti=0.083*(1.935*kappa+1)*Tc;   
            case 3, Kp=0.509*Kc; Td=0.125*Tc; Ti=0.051*(3.302*kappa+1)*Tc;
            case 4
                Kp=(4.437*kappa-1.587)/(8.024*kappa-1.435)*Kc; Ti=0.037*(5.89*kappa+1)*Tc;  Td=0.112*Tc;
        end
    elseif typ == 2
        switch key
            case 2,Kp=(1.892*kappa+0.244)/(3.249*kappa+2.097)*Kc...
                    ; Ti=(0.706*kappa-0.227)/(0.7229*kappa+1.2736)*Tc;   
            case 3, Kp=(4.434*kappa-0.966)/(5.12*kappa+1.734)*Kc...
                    ; Td=0.144*Tc;...
                    Ti=(1.751*kappa-0.612)/(3.776*kappa+1.388)*Tc;
            case 4
                Kp=(4.437*kappa-1.587)/(8.024*kappa-1.435)*Kc; Ti=0.037*(5.89*kappa+1)*Tc;  Td=0.112*Tc;
        end
    end
end
switch key
case 2, nn=Kp*[Ti,1]; dd=[Ti,0];Gc_Sys=tf(nn,dd); H_Sys=[];
case 3
   nn=[Kp*Ti*Td*(N+1)/N, Kp*(Ti+Td/N), Kp];
   dd=Ti*[Td/N,1,0]; Gc_Sys=tf(nn,dd); H_Sys=[];
case 4
   nn=Kp*[Ti,1]; dd=[Ti,0]; dH=Kp*conv([Ti,1],[Td/N,1]);
   nH=[(1+Kp/N)*Ti*Td, Kp*(Ti+Td/N), Kp]; Gc_Sys=tf(nn,dd); H_Sys=tf(nH,dH);
end
