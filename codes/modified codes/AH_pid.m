function [Gc_Sys,H_Sys,Kp,Ti,Td]=AH_pid(key,vars,r)
% Astrom-Hagglund method for PID control design
% r == 1 -> Frequency Response
% r == 0 -> Step Response
% key  = 1 for the PI controller, key = 2 for
% the PID controller  When step response data are available,
% one should specify vars = [K, L, T, N], while vars = [Kc, Tc, Kp, N] are designed for the
% given frequency response data
% Developed by: Sina Masoumi, BSc student, 
% Course: Foundations of Automatic Control Design. Course number: 28255. Fall 1399.
% Department of Mechanical engineering, Sharif University of Technology.
Ti=0; Td=0;
if r == 1
   K=vars(1); Tc=vars(2); Kp=vars(3); N=vars(4); kappa = 1/(K*Kp);
   switch key
   case 1, Kp=0.16*K; Ti = Tc/(1+4.5*kappa);
   case 2, Kp=(0.3-0.1*kappa^4)*K; Ti=0.6/(1+2*kappa)*Tc;...
           Td=0.15*(1-kappa)/(1-0.95*kappa)*Tc; 
   end
else
   K=vars(1); L=vars(2); T=vars(3); N=vars(4); a=K*L/T; 
   switch key
   case 1, Kp=0.15/K + (0.35-L*T/(L+T)^2)/a; Ti=0.35*L + 13*L*T^2/(T^2+12*L*T+7*L^2); 
   case 2, Kp=(0.2+0.45*T/L)/K; Ti=L*(0.4*L+0.8*T)/(L+0.1*T); Td=0.5*L*T/(0.3*L+T); 
   end
end
switch key
case 1, Gc_Sys=tf(Kp*[Ti,1],[Ti,0]); H_Sys=[];
case 2
   dd=Ti*[Td/N,1,0]; nn=[Kp*Ti*Td*(N+1)/N, Kp*(Ti+Td/N), Kp];
   Gc_Sys=tf(nn,dd); H_Sys=[];
end