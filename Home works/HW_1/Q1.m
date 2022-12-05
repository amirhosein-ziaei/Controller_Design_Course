saveVarsMat = load('hw1.mat');

G = saveVarsMat.G; % <1x1 tf> unsupported class

G_delay = saveVarsMat.G_delay; % <1x1 tf> unsupported class

Gd = saveVarsMat.Gd; % <1x1 tf> unsupported class

s = saveVarsMat.s; % <1x1 tf> unsupported class

clear saveVarsMat;
