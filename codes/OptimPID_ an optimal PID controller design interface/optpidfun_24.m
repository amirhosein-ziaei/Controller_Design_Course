function y=optpidfun_24(x)
% OPTPIDFUN_24   An objective function for optimal PID controller tuning
%   It describes the relationship between the PID controller parameters
%   and the chosen performance index.

%   The function is automatically created by PID Controller Optimizer.
%   Date of creation 27-Dec-2021
opt=simset('OutputVariables','y');
assignin('base','Kp',x(1));
assignin('base','Ki',x(2));
assignin('base','Kd',x(3));
try
   [~,~,y_out]=sim('pidctrl_model',[0,200.000000],opt);
catch, y_out=1e10; end
y=y_out(end,1);
