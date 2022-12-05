clc; clear; close all;
s = tf('s');
G_main = 1/(s+1)^8;

[K,L,T]=get_fod(G_main,1);
G = K/(T*s+1);
G.ioDelay = L;

hold on;
legend 'show'
N = 10;
[Kc,pp,wg,wp] = margin(G_main);
Tc = 2*pi/wg;

[Gc1,H1,Kp1,Ti1,Td1] = ziegler_nic(1,[K,L,T,N]);
G_c1 = feedback(Gc1*G_main,1);
G_c1.name = 'P ZN';
step(G_c1);

[Gc2,H2,Kp2,Ti2,Td2] = ziegler_nic(2,[K,L,T,N]);
G_c2 = feedback(Gc2*G_main,1);
G_c2.name = 'PI ZN';
step(G_c2);

[Gc3,H3,Kp3,Ti3,Td3] = ziegler_nic(3,[K,L,T,N]);
G_c3 = feedback(Gc3*G_main,1);
G_c3.name = 'PID ZN';
step(G_c3);

[Gc3,H3,Kp3,Ti3,Td3,beta] = rziegler_nic([K,L,T,N,Kc,Tc]);
step(feedback(Gc3*G_main,H3));