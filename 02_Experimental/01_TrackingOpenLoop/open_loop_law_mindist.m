close all
clear
clc

import casadi.*

load VDP_25s.mat


% define parameter
x  = SX.sym('x',2,1);
u  = SX.sym('u',1,1);
p  = SX.sym('p',2,1);
l  = SX.sym('l',1,1);

nx = length(x);
nu = length(u);

% define time and step size
T = 25;
h = 0.1;

% step-size
N = T/h;
DT = h;

% initial conditions
x0_low = [0.2,0.45]';
x0_up  = x0_low;

% simple bounds
u_low = -1;
u_up  =  1;

x_final = [0.07,0.0]';

R = 0.5;
Q = eye(2);


% dynamics
xdot = [ x(2)
        (1-x(1)^2)*x(2) - x(1) + u];

% control for certain state
u0 = x_final(1)-x_final(2)*(1-x_final(1)^2);


% storage function
V = Function('f',...            % Name
             {p,x},...     % Input variables
             {0.022127*p(1)^2*x(1)^2 + 0.0093277*p(1)^2*x(1)*x(2) + 0.021476*p(1)^2*x(2)^2  ...
              + 2.0362e-08*p(1)*x(1)^3 + 1.1152e-07*p(1)*x(1)^2*x(2) - 8.7296e-08*p(1)*x(1)*x(2)^2  ...
              - 6.1236e-08*p(1)*x(2)^3 + 1.0315*x(1)^4 + 0.34361*x(1)^3*x(2) + 1.6054*x(1)^2 ...
              *x(2)^2 - 0.11699*x(1)*x(2)^3 + 0.87185*x(2)^4 + 0.026562*p(1)*x(1)^2 + 0.076467 ...
              *p(1)*x(1)*x(2) - 0.0010828*p(1)*x(2)^2 - 9.0266e-08*x(1)^3 - 4.7639e-08*x(1)^2 ...
              *x(2) + 1.3461e-07*x(1)*x(2)^2 + 2.2806e-07*x(2)^3 + 0.7776*x(1)^2 + 0.06415 ...
              *x(1)*x(2) + 0.75553*x(2)^2 - 0.25443
});  


h0_low = -inf;
h0_up  = 0;

%% setup optimal control problem

Phi = Function('f',...        
            {x, p},...     
            {V(p,x)}); 


% set up RK45
f = Function('f', {x, u}, {xdot});
X0 = SX.sym('X0', nx,1);
U = SX.sym('U',nu,1);
X = X0;
   
k1 = f(X, U);
k2 = f(X + DT/2 * k1, U);
k3 = f(X + DT/2 * k2, U);
k4 = f(X + DT * k3, U);

X= X + DT/6*(k1 +2*k2 +2*k3 +k4);

F = Function('F', {X0, U}, {X});

% Discretize 
X = casadi.SX.sym('X',nx,2);
U = casadi.SX.sym('U',nu,1);

% vector of decision variables
z = [X(:);U(:)];

% simple bounds on decision variables
lbz =  [-inf(4,1);u_low]; % there are no constraints on the state due to storage function
ubz =  [+inf(4,1);u_up];  % there are no constraints on the state due to storage function

xk   = X(:,1:end-1);
xk1  = X(:,2:end);
uk   = U(:,1:end);

% cost function: stage + end cost
Delta_x = xk1-x_final;
Delta_u = uk-u0;

% J = \alpha * x^T Q x + u^T R *u
J =  p(2)*(Delta_x'*Q*Delta_x) + Delta_u'*Delta_u;

% dynamics as equality constraint
g = xk1 - F(xk,uk);
g = reshape(g,1,size(g,1)*size(g,2));

lbg_dyn = zeros(1,size(X,1)*(size(X,2)-1));
ubg_dyn = lbg_dyn;

% add path constraint here storage function: V(x_(k+1),k*h) 
V = V(p,xk1); % p is varying time step
V = reshape(V,1,size(V,1)*size(V,2));
lbg_cust = h0_low;
ubg_cust = h0_up;

% constraint vecotr
g = [g,V];

% add initial conditions to path constraints 
g   = [xk'  g];

% setup solver
prob   = struct('f', J,...
                'x', z,...
                'g', g,...
                'p',p);


options = struct;                                           
options.print_time = true;

solver = casadi.nlpsol('solver', 'ipopt', prob,options);

 %% open-loop analysis
 alpha_vec = 100; linspace(5,30,20); % alpha from 1 to 100 
 counter   = 0;

 % initialize arrays to store data for analysis
 costFunArr      = zeros(length(alpha_vec),1);
 x_sol_vecArr    = zeros(length(alpha_vec),nx,N+1);
 u_sol_vecArr    = zeros(length(alpha_vec),nu,N);


% analysis for different alpha parameter
for alpha = alpha_vec
    
    costFun = 0;

    % setup solution vector: t = 0 --> k = 0 --> we have N+1 steps because array starts counting at 1
    x_sol_vec = zeros(nx,N+1);
    u_sol_vec = zeros(nu,N);
    
    % first entry of soltuion vector is initial condition
    x_sol_vec(:,1) = x0_low;
   
    % combine path constraints (defect constraints and user defined in-/equality constraints
    lbg = [lbg_dyn,lbg_cust];
    ubg = [ubg_dyn,ubg_cust];

    % first initial points; alternative intial conditions as simple bounds
    lbg = [x0_low' lbg];
    ubg = [x0_low'  ubg];
    
    % initial guess
    z0 = zeros(5,1);
    
    % 1-step optimization open-loop
    for k = 2:N+1

        % solve one-step opt
        [sol]   = solver('x0',  z0,...
                         'p',   [h, alpha],... % dt, alpha (k-1)*h
                         'lbx', lbz,...
                         'ubx', ubz,...
                         'lbg', lbg,...
                         'ubg', ubg);

 
        % extract solution
        z_opt = full(sol.x);
        x_sol = z_opt(1:4);
        u_sol = z_opt(end);

        
        % set initial conditions for next optimization
        lbg(1:2) = x_sol(3:4);
        ubg(1:2) = x_sol(3:4);  
        
        % store states and cost for analysis
        x_sol_vec(:,k)   = x_sol(3:4);
        u_sol_vec(:,k-1)   = u_sol;
        costFun = costFun +  alpha*(x_sol(3:4)-x_final)'*Q*(x_sol(3:4)-x_final) + (u_sol-u0)'*R*(u_sol-u0); 
    
    end

    % store data
    counter                   = counter +1;
    costFunArr(counter)       = costFun;
    x_sol_vecArr(counter,:,:) = x_sol_vec;
    u_sol_vecArr(counter,:,:) = u_sol_vec;
    
    % clear path constraints; are re-initialized in next iteration
    lbg = [];
    ubg = [];

end

%% plotting

% plot sublevel sets and trajectories
pvar x_1 x_2 x_3 t u

domain = [-1 1 -1 1];

% plot level sets of storage function + trajectory
figure('Name','vanderpol-feasible')
clf
hold on
pcontour(subs(Vval, t,0),  0, domain, 'b--')
pcontour(subs(Vval, t,T),  0, domain, 'b-')
plot(x_final(1),x_final(2),'g+')
plot(x0_low(1),x0_low(2),'k+')


% for each alpha one trajectory
[idx_min,~] = find(costFunArr == min(costFunArr));
[idx_max,~] = find(costFunArr == max(costFunArr));

% plot reference trajectory and trajectory minimum and maxmimum cost
for k = [idx_min idx_max]

x_sol_vec = reshape(x_sol_vecArr(k,:,:),nx,N+1);
plot(x_sol_vec(1,:),x_sol_vec(2,:),'-')

end
legend('Sublevel Set at t = 0','Sublevel Set at t = T','Initial Condition', ...
       ['1-step with \alpha = ' num2str(alpha_vec(idx_min))] ,...
       ['1-step with \alpha = ' num2str(alpha_vec(idx_max)) ], 'Location','best')


%% Plot the solution trajectory with best cost function
figure()
tgrid = linspace(0,T,N+1);
clf;
hold on
x_sol_vec = reshape(x_sol_vecArr(idx_min,:,:),nx,N+1);
u_sol_vec = reshape(u_sol_vecArr(idx_min,:,:),nu,N);
plot(tgrid, x_sol_vec(1,:), 'b--')
plot(tgrid, x_sol_vec(2,:), 'r--')
stairs(tgrid, [u_sol_vec nan] , '--k')
legend('x_1','x_2','u')

xlabel('t in seconds')
ylabel('x(t) and u(t)')