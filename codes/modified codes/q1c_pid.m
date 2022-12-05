%1
s=tf('s')
G=(-1.5*s+1)/(s+1)^3;
[K,L,T]=get_fod(G);
N=100;
V=[K,L,T,N];
%ziegler_nicholes
c_z=ziegler_nic(3,V)
hold on
step(feedback(c_z*G,1))
%refined Ziegler-Nichols
% [Kc,p,wc,m]=margin(G)
% [C_rzn,kp,ti,td,beta,H]=rziegler_nic([K L T 10 Kc 2*pi/wc]);
% step(feedback(G*C_rzn,H))
%chr
% 0 overshoot
c_chr_0sh=chr_pid(3,1,[K L T N 0])
step(feedback(c_chr_0sh*G,1))

%20 overshoot
c_chr_20sh=chr_pid(3,1,[K L T N 1])
step(feedback(c_chr_20sh*G,1))

%cohen
c_cohen=cohen_pid(3,V)
step(feedback(c_cohen*G,1))
%WJC
c_wjc=wjcpid(V)
step(feedback(c_wjc*G,1))

%optimum pid
%pid ISTE
c_opt2=opt_pid(3,1,[K L T N 2])
step(feedback(c_opt2*G,1))

legend('ziegler nicholes','CHR 0 overshoot','CHR 20 overshoot','cohen','WJC','opt pid ISTE');
hold off
%pid IST
c_opt1=opt_pid(3,1,[K L T N 1])
%pid IST2E
c_opt3=opt_pid(3,1,[K L T N 3])
%pi_d ISTE
c_opt4=opt_pid(4,1,[K L T N 3])

figure
hold on
step(feedback(c_opt1*G,1))
step(feedback(c_opt2*G,1))
step(feedback(c_opt3*G,1))
step(feedback(c_opt4*G,1))
legend('PID IST','PID ISTE','PID IST2E','PI_D ISTE')