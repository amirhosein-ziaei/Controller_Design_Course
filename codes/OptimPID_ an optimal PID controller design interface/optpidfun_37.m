function y=optpidfun_37(x)
% OPTPIDFUN_37   An objective function for optimal PID controller tuning
%   It describes the relationship between the PID controller parameters
%   and the chosen performance index.

%   The function is automatically created by PID Controller Optimizer.
%   Date of creation 17-Jan-2022
opt=simset('OutputVariables','y');
assignin('base','Kp',x(1));
assignin('base','Ki',x(2));
try
   [~,~,y_out]=sim('pidctrl_model',[0,10.000000],opt);
catch, y_out=1e10; end
y=y_out(end,1);
