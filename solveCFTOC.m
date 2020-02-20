function [feas, xOpt, uOpt, JOpt] = solveCFTOC(nx,nu,Ts, P, Q, R, N, x0, ...
     uL, uU, bf, Af, xDes, uOld)
% This function solves the constraint finite-time optimal control problem
% over the prediction horizon, N. 

% Define state matrix
x = sdpvar(nx,N+1);

% Define input variables
u = sdpvar(nu,N);

% Define objective funciton with summation to N
objective = (xDes - x(:,N+1))'*P*(xDes - x(:,N+1));
for k = 1:N % Matlab indexing needs to be to N instead of N-1
    objective = objective + (x(:,k)-xDes)'*Q*(x(:,k)-xDes)+ u(:,k)'*R*u(:,k);
end

% Implement constraints into yalmip format
% Initialize constraints and initial condition
constraints = [x(:,1)==x0];
for i = 1:N
    constraints = [constraints, uL<=u(:,i)<=uU,...
        x(:,i+1) == modifiedBikeModel(x(:,i), u(:,i), Ts)];
    if i <= N-1       
        constraints = [constraints -0.1<=u(1,i+1)-u(1,i)<=0.1, -deg2rad(8)<=u(2,i+1)-u(2,i)<=deg2rad(8)];
    end
end
constraints = [constraints, -0.1<=u(1,1)-uOld(1)<=0.1, -deg2rad(8)<=u(2,1)-uOld(2)<=deg2rad(8)];

% Set options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','fmincon');

% Terminal set
if isempty(Af)
    % Only penalize for x & y displacement error, not orientation
    constraints = [constraints, x(1:2,N+1)==bf(1:2)];
    Af = 1;
else
    constraints = [constraints, Af*(x(:,N+1)-xDes) <= bf];
end

% Solve and assign solutions to function output
diagnostics = optimize(constraints, objective, options);

% Check feasibility of problem
xOpt = value(x);
uOpt = value(u);
JOpt = value(objective); 

% first 'all' checks columns, second checks all rows
if not( all ( all( uL <= uOpt(:,:) ) ) ) ||  not( all( all( uOpt(:,:) <= uU ) ) )
    feas = 0;
    disp('error 1')
elseif diagnostics.problem ~=0
    disp("ERROR IN YALMIP CFTOC");
    feas = 0;
else
    feas = 1;
end
end