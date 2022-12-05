clc; clear all; close all;
m = 1; g = 10; d = 0.05; L = 1;
s=tf('s');
theta_v = 0.0274/(0.003228*s^2 + 0.003508*s);
r_theta = (m*g*d)/(1.4*L*s^2);
Gp = 0.1*theta_v*r_theta;
Gc =  (s/(s+1))^4 ;
Gc1 = feedback(Gc*Gp,1);
Gc2 = 1/s ; 
G_main = feedback(Gc2 * Gc1 ,1) ;
hold on;
step(G_main,100);

%frequency 
[K,L,T]=get_fod(G_main);
Gc_frequency = K * exp(-L*s)/(T*s+1);

% %opt_app
Gc_optapp  = opt_app(G_main,0,1,1);


figure;
step(G_main , Gc_frequency , Gc_optapp)
legend show;

N = 10;
[Kc,pp,Wg,Wp] = margin ( G_main);
Tc = 2*pi/Wg;

%ziegler_nic
[Gc1,~,~,~,~] = ziegler_nic(3,[K,L,T,N]);
Gc_ZN = 0.5*feedback(Gc1*Gc_frequency,1);

%refined_ziegler_nic
[Gc2,~,~,~,~,~] = rziegler_nic([K,L,T,N,Kc,Tc]);
Gc_RZN = 0.5*feedback(Gc2*G_main,1);

%cohen_coon
[Gc4,~,~,~,~]= cohen_pid (3,1,[K,L,T,N]);
Gc_CC = 0.5*feedback(Gc4*G_main,1);

%cohen_coon_revisited
[Gc5,~,~,~,~]= cohen_pid (3,2,[K,L,T,N]);
Gc_CCR = 0.5*feedback(Gc5*G_main,1);

%AH
[Gc6,~,~,~,~]=astrom_hagglund(2,1,[K,L,T,N]);
Gc_AH =0.5*feedback(Gc6*G_main,1);

%AH_based_frequency
[Gc7,~,~,~,~]=astrom_hagglund(2,2,[K,Kc,Tc,N]);
Gc_AHbf =0.5*feedback(Gc7*G_main,1);

%CHR_0%overshoot
[Gc8,~,~,~,~]=chr_pid(3,1,[K,L,T,N,0]);
Gc_CHR_0OS =0.5*feedback(Gc8*G_main,1);

%CHR_20%overshoot
[Gc9,~,~,~,~]=chr_pid(3,1,[K,L,T,N,1]);
Gc_CHR_20OS =0.5*feedback(Gc9*G_main,1);

%WJC
[Gc10,~,~,~]=wjcpid([K,L,T,N]);
Gc_WJC =0.5*feedback(Gc10*G_main,1);

%optimum_pid
[Gc11,~,~,~,~]=opt_pid(3,1,[K,L,T,N,2]);
Gc_optPID =0.5*feedback(Gc11*G_main,1);

figure;

step(Gc_ZN , Gc_RZN ,Gc_CC , Gc_CCR , Gc_AH , Gc_AHbf , Gc_CHR_0OS , Gc_CHR_20OS , Gc_WJC ,Gc_optPID  )
legend show;