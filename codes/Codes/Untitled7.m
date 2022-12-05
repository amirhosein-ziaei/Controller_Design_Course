clc; clear; close all;
s = tf('s');
G = tf(1,[1 2 0]);

[Kv,L] = get_ipd(G);
G1 = Kv*exp(-L*s)/s;
step(G,G1)
legend 'show'