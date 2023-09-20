close all
clear
clc


%% define variables
x = mpvar('x',2,1);
u = mpvar('u',1,1);

%% dynamics

x_up  = [ 1 1]';
x_low = -[ 1 1]';
u_low = -1;
u_up  = 1;


f = [x(2);
     (1-x(1)^2)*x(2)-x(1)];

gd = [0;1];


f = f+gd*u;


g1 = -(x(1)-x_low(1))*(x_up(1)-x(1));
g2 = -(x(2)-x_low(2))*(x_up(2)-x(2));
hU = [(u-u_low);(u_up-u)];



s1   = polydecvar('s1' ,monomials(x,0:2));
s2   = polydecvar('s2' ,monomials(x,0:2));

s3   = polydecvar('s3' ,monomials(x,0:2));
s4   = polydecvar('s4' ,monomials(x,0:2));

s31 = polydecvar('c31' ,monomials(x,0:2));
s32  = polydecvar('c31' ,monomials(x,0:2));

s41 = polydecvar('c41' ,monomials(x,0:2));
s42 = polydecvar('c42' ,monomials(x,0:2));
 
V = polydecvar('V',monomials(x(1),1:1));
k = polydecvar('k',monomials(x,1));

[A,B] = plinearize(f, x, u);
[K,P] = lqr(A, B, eye(2), 2.5);

Pinit = lyap(A-B*K,eye(2)); % P calculated for scaled dynamics
% l     = x'*P*x-1;           % here x is already scaled
Vval  = x'*x;           % here x is already scaled

sopt = sosoptions;
sopt.solver = 'mosek';
% tic
% for iter = 1:maxIter
 
    % solve for multiplier and control law
    sosconk = [s1 >= 0;
               s2 >= 0; 
               s31 >= 0; 
               s32 >= 0; 
               s41 >= 0; 
               s42 >= 0;
               s1*Vval - jacobian(Vval,x)*subs(f,u,k) >= 0;
               s2*Vval + jacobian(Vval,x)*subs(f,u,k) >= 0;
               s3*Vval + subs(f,u,k) >= 0;
               s31*Vval - g1 >= 0;
               s32*Vval - g2 >= 0;
               s41*Vval + hU(1) >= 0;
               s42*Vval + hU(2) >= 0;];

   [info,dopt] = sosopt(sosconk,x,sopt)


   s1val = subs(s1,dopt);
   s2val = subs(s2,dopt);
   s31val = subs(s31,dopt);
   s32val = subs(s32,dopt);
   s41val = subs(s41,dopt);
   s42val = subs(s42,dopt);





% end