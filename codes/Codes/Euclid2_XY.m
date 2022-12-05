function [Ns, Ms, X, Y] = Euclid2_XY (G)
% input: transfer function G
% Output N(s), M(s), X(s), Y(s): G(s)=N/M,  NX+MY=1.
% la=lambda , here: lambda=1/(s+alpha), % you can change alpha in line 8.
% you may need to use -1*X, -1*Y
% course: Foundamentals of Automatic Control Design. 28-255 , Term: 1398-99-2
% Sharif University of Technology, Department of Mechanical Engineering, Tehran, Iran.
% Prepared by: Mohammad Jandaghi Anaraki, S.N.: 95104555
% 1399/03/30
% The file has been modified by Sina Masoumi, BSc student, Fall 1399.
%
[~,k] = size(G.Denominator{1});
s = tf('s');
k = k-1;
alpha = 1;
numerator = tf(G.Numerator{1},1);
denominator = tf(G.Denominator{1},1);
Ns = numerator/(s+alpha)^k;
Ms = denominator/(s+alpha)^k;
syms x landa
Nx = simplify(poly2sym(Ns.Numerator{1})/poly2sym(Ns.Denominator{1}, x));
Mx = simplify(poly2sym(Ms.Numerator{1})/poly2sym(Ms.Denominator{1}, x));
x = 1/landa - alpha;
nl = simplify(subs(Nx));
ml = simplify(subs(Mx));
g = simplify(nl/ml);
[nl, ml] = numden(g);
n = sym2poly(nl);
m = sym2poly(ml);
sn = size(n);
sm = size(m);
if sn(2) < sm(2)
    f = n;
    n = m;
    m = f;
    f = Nx;
    Nx = Mx;
    Mx = f;
end
[q1, r1]=deconv(n,m);
r_new = nonzeros(r1)';
b = m;
j = 0;
Q_matrix = zeros(length(m)-1,max(2,length(n)-length(m)+1));
R_matrix = zeros(length(m)-1,length(m)-1);

syms s 
lambda = 1/(s+alpha);

if (length(n) == length (m))
    Q_matrix(1,end) = q1;
    R_matrix(1,:) = r_new;
else
    Q_matrix(1,:) = q1;
    R_matrix(1,:) = r_new;
end

while (length(r_new)>1)
    c = r_new;
    [q_new, r_new] = deconv(b,r_new);
    b = c;
    r_new = nonzeros(r_new)' ;
    Q_matrix(j+2,end-1:end) = q_new ;
    R_matrix(j+2,j+2:end) = r_new ;
    j = j+1 ;
end


if (j == 0)
    x = 1/R_matrix(end);
    y = -Q_matrix(1,:)/R_matrix(end);
    
    X = 0;
    Y = 0;
    for i = 1:length(x)
        X = X + x(i)*lambda^(length(x)-i);
    end
    
    for i = 1:length(y)
        Y = Y + y(i)*lambda^(length(y)-i);
    end
    
    X = vpa(simplify(expand(X)),5);
    Y = vpa(simplify(expand(Y)),5);
    
end

if(j == 1)
    one_q1q2 = zeros(1,length(Q_matrix(1,:))+length(Q_matrix(2,:))-1);
    one_q1q2(end) = 1 ;
    
    x = -Q_matrix(2,:)/R_matrix(end) ;
    y = (conv(Q_matrix(1,:),Q_matrix(2,:))+ one_q1q2) / R_matrix(end) ;

    X = 0;
    Y = 0;
    for i = 1:length(x)
        X = X + x(i)*lambda^(length(x)-i);
    end
    
    for i = 1:length(y)
        Y = Y + y(i)*lambda^(length(y)-i);
    end
    
    X = vpa(simplify(expand(X)),5);
    Y = vpa(simplify(expand(Y)),5);
    
end
    
    if(j == 2)
        
        one_q2q3 = zeros(1,length(Q_matrix(2,:))+length(Q_matrix(3,:))-1);
        one_q2q3(end) = 1;
        
        x = (one_q2q3 +conv(Q_matrix(2,:),Q_matrix(3,:)))/R_matrix(end);
        
        a = (one_q2q3 +conv(Q_matrix(2,:),Q_matrix(3,:)));
        b = conv(Q_matrix(1,:),a);
        q3 = [zeros(1,length(b)-2),Q_matrix(3,:)];
        y = (-q3 - b) / R_matrix(end);
        X = 0;
        Y = 0;
        for i = 1:length(x)
            X = X + x(i)*lambda^(length(x)-i);
        end
        
        for i = 1:length(y)
        Y = Y + y(i)*lambda^(length(y)-i);
        end
        
        X = vpa(simplify(expand(X)),5);
        Y = vpa(simplify(expand(Y)),5);
        
    end
syms s
x = s;
Ns = subs(Nx);
Ms = subs(Mx);
vpa(simplify(Ns*X+Ms*Y),2)
% if vpa(simplify(Ns*X+Ms*Y)) < 0
%     X = -1*X;
%     Y = -1*Y;
end
