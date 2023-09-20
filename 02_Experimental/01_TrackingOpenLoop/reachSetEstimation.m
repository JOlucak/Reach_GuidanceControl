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


% [Dx,rx,Du,ru] = calcScales(x_up,x_low,u_up,u_low);
% 
Dx = eye(2);
Du = 1;

rx = zeros(2,1);
ru = 0;

xs = Dx*x+rx;
xu = Dx^-1*(x-rx);

us = Du*u+ru;
uu = Du^-1*(u-ru);

f = [x(2);
     (1-x(1)^2)*x(2)-x(1)];

gd = [0;1];


x_ups  = Dx*x_up+rx;
x_lows = Dx*x_low+rx;

u_ups  = Du*u_up+ru;
u_lows = Du*u_low+ru;


%% Terminal Penalty & Constraints
% terminal LQR
[A,B] = plinearize(subs(Dx*(f+gd*u),[x;u],[xu;uu]), x, u);
[K,P] = lqr(A, B, eye(2), 2.5);

Pinit = lyap(A-B*K,eye(2)); % P calculated for scaled dynamics
l     = x'*P*x-1;           % here x is already scaled
Vval  = x'*P*x-1;           % here x is already scaled

g1 = -(x(1)-x_low(1))*(x_up(1)-x(1));
g2 = -(x(2)-x_low(2))*(x_up(2)-x(2));


hU = [(u-u_low);(u_up-u)];
%% Problem Definitions
T     = 25;
t1    = 0;
t2    = T;

deg_V = 2;

maxIter = 10;


%% define SOS problem
pvar t

hT = (t-t1)*(t2-t);


s1  = sosdecvar('c1',monomials([t;x],0:4));
s2  = sosdecvar('c21',monomials([t;x],0:2));

s31 = sosdecvar('c31',monomials([t;x],0:2));
s32 = sosdecvar('c32',monomials([t;x],0:2));
s41 = sosdecvar('c41',monomials([t;x],0:2));
s42 = sosdecvar('c42',monomials([t;x],0:2));

s51  = sosdecvar('c51',monomials([t;x],0:2));
s52  = sosdecvar('c52',monomials([t;x],0:2));
s61  = sosdecvar('c61',monomials([t;x],0:2));
s62  = sosdecvar('c62',monomials([t;x],0:2));
s7   = sosdecvar('c7' ,monomials(x,0:2));
s8   = sosdecvar('c8' ,monomials(x,0:2));
 
V = polydecvar('V',monomials([t;x],0:5));
k = polydecvar('k',monomials([t;x],0:2));



domain = [-1 1 -1 1];
pcontour(l,0,domain,'r')
hold on
line([x_lows(1) ,x_ups(1)],[x_ups(2) ,x_ups(2)],'Color','black','LineStyle','--')
line([x_lows(1) ,x_ups(1)],[x_lows(2) ,x_lows(2)],'Color','black','LineStyle','--')
line([x_lows(1) ,x_lows(1)],[x_lows(2) ,x_ups(2)],'Color','black','LineStyle','--')
line([x_ups(1) ,x_ups(1)],[x_lows(2) ,x_ups(2)],'Color','black','LineStyle','--')
hold on
pcontour(Vval,0,domain,'b--')

sopt = sosoptions;
% sopt.solver = 'mosek';
tic
for iter = 1:maxIter
 
    % solve for multiplier and control law
    sosconk = [s1 >= 0;
               s2 >= 0; 
               s31 >= 0; 
               s32 >= 0; 
               s41 >= 0; 
               s42 >= 0;
               s51 >= 0; 
               s52 >= 0; 
               s61 >= 0; 
               s62 >= 0;...
               subs(hU(1),u,k) + s51*Vval - s61*hT >= 0;...
               subs(hU(2),u,k) + s52*Vval - s62*hT >= 0;...
               s1*Vval - s2*hT - jacobian(Vval,t) - jacobian(Vval,x)*f - jacobian(Vval,x)*gd*k >= 0;...
               s31*Vval - s41*hT - g1 >= 0;
               s32*Vval - s42*hT - g2 >= 0];

[info,dopt] = sosopt(sosconk,[x;t],sopt);

if ~info.feas
    disp(['h step infeasible in iteration: ' num2str(iter)]); 
    break
else
    disp(['h step feasible in iteration:' num2str(iter)])
    s1val  = subs(s1,dopt);
    s2val  = subs(s2,dopt);

    s31val  = subs(s31,dopt);
    s32val  = subs(s32,dopt);

    s41val  = subs(s41,dopt);
    s42val  = subs(s42,dopt);

    s51val = subs(s51,dopt);
    s52val = subs(s52,dopt);
    s61val = subs(s61,dopt);
    s62val = subs(s62,dopt);
    kval   = subs(k,dopt);

end

 % solve for storage function
sosconV = [s7 >= 0;...
           s1val*V - s2val*hT - jacobian(V,t) - jacobian(V,x)*f - jacobian(V,x)*gd*kval >= 0;...
           s31val*V - s41val*hT - g1 >= 0;
           s32val*V - s42val*hT - g2 >= 0;...
           subs(hU(1),u,kval) + s51val*V - s61val*hT >= 0;...
           subs(hU(2),u,kval) + s52val*V - s62val*hT >= 0;...
           subs(V,t,t2) - l >= 0;...
           s7*subs(Vval,t,t1)-subs(V,t,t1) >= 0];

sopt = sosoptions;
[info,dopt] = sosopt(sosconV,[x;t],sopt);
   
if ~info.feas
    disp(['V step infeasible in iteration:' num2str(iter)]); 
    break
else
    disp(['V step feasible in iteration:' num2str(iter)])
    Vval = subs(V,dopt);
end
end
toc

Vval
% Vu = subs(Vval,x,xu)
% ku = subs(kval,x,xu)
%% plotting
domain = [-1 1 -1 1];
figure()
pcontour(subs(Vval,t,0),0,domain,'b--')
hold on
pcontour(subs(Vval,t,T),0,domain,'b')
hold on
pcontour(subs(l,x,x),0,domain,'r')
hold on
hold on
line([x_low(1) ,x_up(1)],[x_up(2) ,x_up(2)],'Color','black','LineStyle','--')
line([x_low(1) ,x_up(1)],[x_low(2) ,x_low(2)],'Color','black','LineStyle','--')
line([x_low(1) ,x_low(1)],[x_low(2) ,x_up(2)],'Color','black','LineStyle','--')
line([x_up(1) ,x_up(1)],[x_low(2) ,x_up(2)],'Color','black','LineStyle','--')
legend('BRS','Inner Terminal','Terminal Set','State Constraint')

% x_up  = Dx^-1*(x_up-rx);
% x_low = Dx^-1*(x_low-rx);
% 
% u_up  = Du^-1*(u_up-ru);
% u_low = Du^-1*(u_low-ru);
% 
% figure()
% Vvalu = subs(Vval,x,inv(Dx)*(x-rx));
% pcontour(subs(Vvalu,t,0),0,domain,'b--')
% hold on
% pcontour(subs(Vvalu,t,T),0,domain,'b')
% hold on
% pcontour(subs(l,x,(Dx)*x+rx),0,domain,'r')
% hold on
% hold on
% line([x_low(1) ,x_up(1)],[x_up(2) ,x_up(2)],'Color','black','LineStyle','--')
% line([x_low(1) ,x_up(1)],[x_low(2) ,x_low(2)],'Color','black','LineStyle','--')
% line([x_low(1) ,x_low(1)],[x_low(2) ,x_up(2)],'Color','black','LineStyle','--')
% line([x_up(1) ,x_up(1)],[x_low(2) ,x_up(2)],'Color','black','LineStyle','--')
