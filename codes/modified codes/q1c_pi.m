%1
s=tf('s')
G=(-1.5*s+1)/(s+1)^3;
[K,L,T]=get_fod(G);
N=100;
V=[K,L,T,N];
%ziegler_nicholes
c_z=ziegler_nic(2,V)
hold on
step(feedback(c_z*G,1))
%refined Ziegler-Nichols
% [Kc,p,wc,m]=margin(G)
% [C_rzn,kp,ti,td,beta,H]=rziegler_nic([K L T 10 Kc 2*pi/wc]);
% step(feedback(G*C_rzn,H))
%chr
% 0 overshoot
c_chr_0sh=chr_pid(2,1,[K L T N 0])
step(feedback(c_chr_0sh*G,1))

%20 overshoot
c_chr_20sh=chr_pid(2,1,[K L T N 1])
step(feedback(c_chr_20sh*G,1))


%WJC
c_wjc=wjcpid(V)
step(feedback(c_wjc*G,1))

%optimum pid
%pid ISTE
legend('ziegler nicholes','CHR 0 overshoot','CHR 20 overshoot','WJC');
hold off
