clc; clear; close all;
s = tf('s');
G1 = 1/(s+1)^8;

[K,L,T]=get_fod(G1);
G1_freq = K*exp(-L*s)/(T*s+1);

[K,L,T]=get_fod(G1,1);
G1_tf = K*exp(-L*s)/(T*s+1);
step(G1,G1_freq,G1_tf)
legend 'show'
