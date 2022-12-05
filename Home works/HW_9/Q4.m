clc;clear all;
s=tf('s');
G = 0.1*s/((s-1)*(s+0.1));
[Ns, Ms, X, Y] = Euclid2_XY (G);
