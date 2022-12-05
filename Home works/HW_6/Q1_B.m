clc; clear all ; close all;
s = tf('s');
G_main = (-s+3)/((s+1)*(s+2)*(s^2 + 2*s + 4));

%transfer_function
[K,L,T] = get_fod(G_main,1);
G_tf = K*exp(-L*s)/(T*s + 1);

[Kc,pp,Wg,Wp] = margin ( G_main);
Tc = 2*pi/Wg;

N = 10;
rb=0.5;
pb=0.5;
legend ' show ' ;
hold on;

% %ziegler_nic
[Gc1,~,~,~,~] = ziegler_nic(3,[K,L,T,N]);
Gc_ZN = feedback(Gc1*G_main,1);


% %refined_ziegler_nic
[Gc2,~,~,~,~,~] = rziegler_nic([K,L,T,N,Kc,Tc]);
Gc_RZN = feedback(Gc2*G_main,1);


%modified_ziegler_nic
[Gc3,~,~,~,~] = ziegler_nic(3,[K,Tc,rb,pb,N]);
Gc_modified_ZN = feedback(Gc3*G_main,1);


%cohen_coon
[Gc4,~,~,~,~]= cohen_pid (3,1,[K,L,T,N]);
Gc_CC = feedback(Gc4*G_main,1);


%cohen_coon_revisited
[Gc5,~,~,~,~]= cohen_pid (3,2,[K,L,T,N]);
Gc_CCR = feedback(Gc5*G_main,1);


%AH
[Gc6,~,~,~,~]=astrom_hagglund(2,1,[K,L,T,N]);
Gc_AH =feedback(Gc6*G_main,1);


%AH_based_frequency
[Gc7,~,~,~,~]=astrom_hagglund(2,2,[K,Kc,Tc,N]);
Gc_AHbf =feedback(Gc7*G_main,1);


%CHR_0%overshoot
[Gc8,~,~,~,~]=chr_pid(3,1,[K,L,T,N,0]);
Gc_CHR_0OS =feedback(Gc8*G_main,1);


%CHR_20%overshoot
[Gc9,~,~,~,~]=chr_pid(3,1,[K,L,T,N,1]);
Gc_CHR_20OS =feedback(Gc9*G_main,1);


%WJC
[Gc10,~,~,~]=wjcpid([K,L,T,N]);
Gc_WJC =feedback(Gc10*G_main,1);

%optimum_pid
[Gc11,~,~,~,~]=opt_pid(3,2,[K,L,T,N,2]);
Gc_optPID =feedback(Gc11*G_main,1);


%optimum_pi_d
[Gc12,H_PI_D,~,~,~]=opt_pid(4,2,[K,L,T,N,2]);
Gc_optPI_D =feedback(Gc12*G_main,H_PI_D);


step (Gc_ZN ,Gc_RZN,Gc_modified_ZN,Gc_CC,Gc_CCR,Gc_AH,Gc_AHbf,Gc_CHR_0OS,Gc_CHR_20OS,Gc_WJC,Gc_optPID,Gc_optPI_D);
