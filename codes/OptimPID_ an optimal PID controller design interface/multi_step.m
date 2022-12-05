function [sys,x0,str,ts]=multi_step(t,x,u,flag,tTime,yStep)
switch flag,
   case 0 % 调用初始化过程
       [sys,x0,str,ts] = mdlInitializeSizes;
   case 3  % 计算输出信号，生成多阶跃信号
       sys = mdlOutputs(t,tTime,yStep);
   case {1, 2, 4, 9} % 未使用的 flag 值
       sys = [];
   otherwise % 错误信息处理
       error(['Unhandled flag = ',num2str(flag)]);
end;
% when flag=0 时，进行初始化处理
function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes; % 调入初始化的模版
sizes.NumContStates = 0; sizes.NumDiscStates = 0; % 无连续、离散状态
sizes.NumOutputs = 1; sizes.NumInputs = 0; % 系统的输入和输出路数
sizes.DirFeedthrough = 0;% 输入信号不直接传输到输出
sizes.NumSampleTimes = 1;% 单个采样周期
sys = simsizes(sizes); % 初始化
x0 = []; str = []; % 系统的初始状态为空向量
ts = [0 0]; % 假设模块为连续模块
%  flag=3 时，计算输出信号
function sys = mdlOutputs(t,tTime,yStep)
i=find(tTime<=t); sys=yStep(i(end));
