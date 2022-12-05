clc; clear all; close all;
m = 1; g = 10; d = 0.05; L = 1;
s=tf('s');
theta_v = 0.0274/(0.003228*s^2 + 0.003508*s);
r_theta = (m*g*d)/(1.4*L*s^2);
Gp = 0.1*theta_v*r_theta;
step(feedback(Gp,1))
figure;
rlocus(Gp)