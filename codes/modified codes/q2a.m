s=tf('s');
G=12*(s^2-3*s+6)/((s+1)*(s+5)*(s^2+3*s+6)*(s^2+s+2))
[K,L,T]=get_fod(G);
fotd=K/(T*s+1)*exp(-s*L);
OP=opt_app(G,0,1,1)
step(G,fotd,OP)
legend('G','getfod','opt_app')