clc; clear all; close all; 
s = tf('s');
G_main = (-s+3)/((s+1)*(s+2)*(s^2 + 2*s + 4));

%frequency
[K,L,T] = get_fod(G_main);
G_frequency = K*exp(-L*s)/(T*s + 1);

%transfer_function
[K1,L1,T1] = get_fod(G_main,1);
G_tf = K1*exp(-L1*s)/(T1*s + 1);

%opt_app
G_opt_app  = opt_app(G_main,0,1,1);

legend ' show ' 
step ( G_main , G_frequency , G_tf , G_opt_app );