%% Cassie Hovershoe MPC for Trajectory Tracking
close all 
clear all
clc

%% Initalize Environment
% add paths
startup;

% Load Cassie model 
if exist('./mat/cassie_hovershoe_model.mat')
    load('./mat/cassie_hovershoe_model.mat') ;
else
    model = create_cassie_featherstone_model() ;
end

% Initial configuration: 62-by-1
x0_full = getInitialStateCassieHovershoe(model);

% Get control parameters
params = hovershoeParams(model);

% Define timestep
externalForce_fun = @ExternalForce;

%% Model Predicitive Control
% Hyper Parameters
M = 175; N = 30; Ts = 1.0; 

nx = 3; % number of states
nu = 2; % number of inputs

Q = [10  0  0;              % State Penalty
      0 10  0;
      0  0  0];
P = 1*eye(nx,nx);           % Terminal State Penalty
R = eye(nu,nu);             % Input Penality


% Create Sinusoidal Reference Track
x_coor = linspace(0,12*pi,20);
y_coor = 1*cos(0.25*x_coor) - 1;
plot(x_coor,y_coor,'o');

xDes_all(1,:) = x_coor;
xDes_all(2,:) = y_coor;
xDes_all(3,:) = 0;
xDes = xDes_all(:,2);
counter = 1;


% Get position of center of mass
r0_com = compute_COM_pos(model, x0_full(1:model.NB))' ;

x0 = [r0_com(1); % initial X position of the COM
      r0_com(2); % initial Y position of the COM
      0];

% Input Box Constraints
uL = [0, deg2rad(-24)]';    % Umin
uU = [0.2, deg2rad(24)]';   % Umax

% Terminal set constraints
Af = [1,  0,  0;
     -1,  0,  0;
      0,  1,  0;
      0, -1,  0];
bf = 0.2*[1; 1; 1; 1];

% ODE Options
time_inter = [0 Ts]; 
odeopts = odeset('Events', @falldetect_MPC,'AbsTol',1e-3,'RelTol',1e-6);
externalForce_fun = @ExternalForce;

% Init storing vectors
xOpt_com(:,1) = x0; % Reduced states for simplfied model for MPC
xOpt_full(1,:) = x0_full'; % Full 62 states for simulation
xPred = zeros(nx,N+1,M);
predErr = zeros(nx,M-N);
feas = false([1,M]);
t_vec_full(1) = 0;
t_vec_store = [];
x_vec_store = [];
r0_com_store = [];
savevec = [];
uOld = [0; 0];

for t = 1:M
    fprintf('Solving simstep: %i\n',t)
    
    [feas(t), x, u] = solveCFTOC(...
        nx,nu,Ts,P,Q,R,N,xOpt_com(:,t),uL,uU,bf,Af,xDes,uOld);
    
    if ~feas(t)
        return;
        disp("Error in CFTOC");
    end
    
    % Save closed loop trajectory and first control input for simulation
    uOpt(:,t) = u(:,1); 
    uStar = u(:,1); % which is first element of each CFTOC solution 
    
    % Apply uStar to the system and simulate 1 step (Ts)
    [t_vec, x_vec, u_theta_vec, u_psi_vec] = ...
        ode15s(@cassieHovershoeEOM_MPC, time_inter, xOpt_full(t,:)',...
        odeopts, model, params, externalForce_fun, uStar, N,Ts, uOld);
    
    % Get new states, set as the inital states, and repeat
    uOld = uStar;
    t_vec_full(t+1) = t*Ts;
    xOpt_full(t+1,:) = x_vec(end,:);
    
    % Recompute where the COM has traveled
    r0_com_new = compute_COM_pos(model, x_vec(end,1:model.NB))'; 
    
    xOpt_com(1:2,t+1) = r0_com_new(1:2);    % 1:2 are X and Y positon of COM
    xOpt_com(3,t+1) = x_vec(end,4);         % 4th State is yaw
    
    t_vec = t_vec + (t-1)*Ts;
    t_vec_store = [t_vec_store; t_vec];
    x_vec_store = [x_vec_store; x_vec];
    
    % check if norm(xDes - COMposition) at every step is close enough to
    % move to the next WP
    checknorm = [];
    for j = 1:size(x_vec,1)
        compos = compute_COM_pos(model, x_vec(j,1:model.NB));
        checknorm(j) = norm(xDes(1:2) - compos(1:2));
    end
    
    savevec = [savevec; checknorm'];
    % Transition to next WP if current position and current WP are within
    % 0.5m
    if min(checknorm) < 0.5
        counter = counter + 1;
        xDes = xDes_all(:,counter);
    end
end

% display the optimal path
disp(xOpt_com);
disp(uOpt);

%% Visualize
% Animation
% 3D model animation:
stateData = getVisualizerState_CassieHover(x_vec_store, model);
vis = CassieHoverVisualizer(t_vec_store, stateData);

%% Plots
% Plot Hovershoes & COM trajectory
figure ;
    plot(xOpt_full(1:177,49),xOpt_full(1:177,50),'LineWidth',3); hold on;
    plot(xOpt_full(1:177,56),xOpt_full(1:177,57),'LineWidth',3); hold on;
    plot(xOpt_com(1,1:177), xOpt_com(2,1:177),'LineWidth',3);
    xlabel('X-axis (m)','FontSize',30) ; grid ;ylabel('Y-axis (m)','FontSize',30);
    title('Turning Trajectory','FontSize',30) ;
    title('x-y plot','FontSize',30);
    
    %axis equal
    
% plot waypoints
for i = 1:13
    hold on
    plot(xDes_all(1,i), xDes_all(2,i),'k*','MarkerSize',12,'lineWidth',2)
end
legend('Left Hovershoe','Right Hovershoe','Cassie COM','Wayponts','FontSize',30,'Location','NorthEastOutside');

% Plot Cassie dyaw    
figure ;
plot(t_vec_full, xOpt_full(:,24));
xlabel('Time (s)') ; grid ;ylabel('dYaw (rad/s)') ; title('Yaw Derivative') ;
title('dyaw');

% Plot Cassie velocity    
figure ;
plot(t_vec_full(:,:), (xOpt_full(:,21).^2 + xOpt_full(:,22).^2).^0.5);
xlabel('Time (s)'); grid; ylabel('Velocity (m/s)'); title('Speed') ;
title('Hovershoe Velocity');