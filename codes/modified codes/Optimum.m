function [Gc, Kp, Ti, Td, H] = Optimum(key, typ, vars)
%
% Design PID controller using the optimal setting algorithm.
% This function supposed to be equal to function opt_pid.m.
% The disturbance rejection part has been added to the main file.
% This function has been prepared by:
% Erfan Etesami, BSc student 96106214, Dept of Mech Eng, Sharif U of T, Dec 2020, Dey 1399.
% 
% key == 2 -> PI, key == 3 -> PID, key == 4 -> PI-D
% typ == 1 -> set point, typ == 2 -> disturbance rejection
% length(vars) == 5 -> [K, L, T, N, iC]
% length(vars) == 7 -> [K, L, T, N, Kc, Tc, kappa]
% iC == 1 -> ISE, iC == 2 -> ISTE, iC == 3 -> IST2E

Ti = [];
Td = [];
H = 1;

K = vars(1); 
L = vars(2);
T = vars(3);
N = vars(4);

if length(vars) == 5 
    iC = vars(5);
    
    if typ == 1 % set point
        switch key
            case 2   % PI
                % Table 6.5
                PIDtab = [0.980, 0.712, 0.569, 1.072, 0.786, 0.628;
                          -0.892, -0.921, -0.951, -0.560, -0.559, -0.583;
                          0.690, 0.968, 1.023, 0.648, 0.883, 1.007;
                          -0.155, -0.247, -0.179, -0.114, -0.158, -0.167];
            case 3  % PID
                % Table 6.6
                PIDtab = [1.048, 1.042, 0.968, 1.154, 1.142, 1.061;
                          -0.897, -0.897, -0.904, -0.567, -0.579, -0.583;
                          1.195, 0.987, 0.977, 1.047, 0.919, 0.892;
                          -0.368, -0.238, -0.253, -0.220, -0.172, -0.165;
                          0.489, 0.385, 0.316, 0.490, 0.384, 0.315;
                          0.888, 0.906, 0.892, 0.708, 0.839, 0.832];
            case 4  % PI-D
                % Table 6.7
                PIDtab = [1.260, 1.053, 0.942, 1.295, 1.120, 1.001;
                          -0.887, -0.930, -0.933, -0.619, -0.625, -0.624;
                          0.701, 0.736, 0.770, 0.661, 0.720, 0.754;
                          -0.147, -0.126, -0.130, -0.110, -0.114, -0.116;
                          0.375, 0.349, 0.308, 0.378, 0.350, 0.308;
                          0.886, 0.907, 0.897, 0.756, 0.811, 0.813];
        end
        % end switch
        
        ii = 0;
        tt = L/T;
        if tt>1
            ii = 3;
        end
        
        a1 = PIDtab(1, ii+iC);
        b1 = PIDtab(2, ii+iC);
        a2 = PIDtab(3, ii+iC);
        b2 = PIDtab(4, ii+iC);
        
        Kp = (a1/K) * (tt^b1);
        Ti = T / (a2 + b2*tt);
        
        if key==3 || key==4
            a3 = PIDtab(5, ii+iC);
            b3 = PIDtab(6, ii+iC);
            Td = a3 * T * (tt^b3);
        end
        
    elseif typ == 2 % disturbance rejection 
        switch key
            case 2   % PI
                % Table 6.8
                PIDtab = [1.279, 1.015, 1.021, 1.346, 1.065, 1.076;
                          -0.945, -0.957, -0.953, -0.675, -0.673, -0.648;
                          0.535, 0.667, 0.629, 0.552, 0.687, 0.650;
                          0.586, 0.552, 0.546, 0.438, 0.427, 0.442];
            case 3  % PID
                % Table 6.9
                PIDtab = [1.473, 1.468, 1.531, 1.524, 1.515, 1.592;
                          -0.970, -0.970, -0.960, -0.735, -0.730, -0.705;
                          1.115, 0.942, 0.971, 1.130, 0.957, 0.957;
                          0.753, 0.725, 0.746, 0.641, 0.598, 0.597;
                          0.550, 0.443, 0.413, 0.552, 0.444, 0.414;
                          0.948, 0.939, 0.933, 0.851, 0.847, 0.850];
        end
        % end switch
        
        ii = 0;
        tt = L/T;
        if tt>1
            ii = 3;
        end
        
        a1 = PIDtab(1, ii+iC);
        b1 = PIDtab(2, ii+iC);
        a2 = PIDtab(3, ii+iC);
        b2 = PIDtab(4, ii+iC);
        
        Kp = (a1/T) * (tt^b1);
        Ti = (T/a2) * (tt^b2);
        
        if key == 3
            a3 = PIDtab(5, ii+iC);
            b3 = PIDtab(6, ii+iC);
            Td = a3 * T * (tt^b3);
        end
    end
    
elseif length(vars) == 7 
    Kc = vars(5); 
    Tc = vars(6); 
    kappa = vars(7); 
    
    if typ == 1 % set point
        switch key
            case 2  % PI
                Kp = Kc * (4.264-0.148*kappa) / (12.119-0.432*kappa);
                Ti = 0.083 * Tc * (1.935*kappa+1);
            case 3  % PID
                Kp = 0.509 * Kc;
                Ti = 0.051 * Tc * (3.302*kappa+1);
                Td = 0.125 * Tc;
            case 4  % PI-D
                Kp = Kc * (4.437*kappa-1.587) / (8.024*kappa-1.435);
                Ti = 0.037 * Tc * (5.89*kappa+1);
                Td = 0.112 * Tc;
        end
    elseif typ == 2 % disturbance rejection
        switch key
            case 2  % PI
                Kp = Kc * (1.892*kappa+0.244) / (3.249*kappa+2.097);
                Ti = Tc * (0.706*kappa-0.227) / (0.7229*kappa+1.2736);
            case 3  % PID
                Kp = Kc * (4.434*kappa-0.966) / (5.12*kappa+1.734);
                Ti = Tc * (1.751*kappa-0.612) / (3.776*kappa+1.388);
                Td = 0.144 * Tc;
        end
    end
end

% PID Controller Design
switch key
    case 2  % PI
%         nn = Kp * [Ti, 1];
%         dd = [Ti, 0];
%         Gc = tf(nn, dd);
        s = tf('s');
        Gc = Kp * (1 + 1/Ti/s);
    case 3  % PID
%         nn = Kp * [Ti*Td*(N+1)/N, (Ti+Td/N), 1];
%         dd = Ti * [Td/N, 1, 0];
%         Gc = tf(nn, dd);
        s = tf('s');
        Gc = Kp * (1 + (1/Ti/s) + (Td*s)/(1+Td*s/N));
    case 4  % PI-D
%         nn = Kp * [Ti, 1];
%         dd = [Ti, 0];
%         Gc = tf(nn, dd);
        s = tf('s');
        Gc = Kp * (1 + 1/Ti/s);
%         nH = [(1+Kp/N)*Ti*Td, Kp*(Ti+Td/N), Kp];
%         dH = Kp * conv([Ti, 1], [Td/N, 1]);
%         H = tf(nH, dH);
        Gd = (Td*s) / (1+Td*s/N);
        H = 1 + Gd / Gc;
end

end