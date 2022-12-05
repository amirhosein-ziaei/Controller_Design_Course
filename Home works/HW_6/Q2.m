clc; clear all ; close all;
s = tf('s');
G_main = (-s+3)/(s*(s+2)*(s^2 + 2*s + 4));
N=10;
G1 = G_main*s;
Gc= opt_app(G1,0,1,1);
G_opt_app=Gc/s;

[Kv,L]= get_ipd(G_main);
G_ipd=tf(Kv,[1,0]);
G_ipd.iodelay = L;

step ( G_main,G_ipd,G_opt_app);
G_c = feedback(G_main,1);

[Gc_PD,~,~,~] = ipdtctrl(1,2,Kv,L,N);
G_PD = feedback(Gc_PD*G_main,1);

[Gc_PID,~,~,~] = ipdtctrl(2,2,Kv,L,N);
G_PID = feedback(Gc_PID*G_main,1);

step ( G_c,G_PD,G_PID);

