clc; clear; close all;

%% =================== Parameter Settings ===================
m_d = 1;          % Admittance mass
b_d = 6;          % Admittance damping
m = 0.5;          % Actual system mass

mu = 0.8;         % Coulomb friction coefficient
b = 0.4;          % Viscous friction coefficient

Ad = 10;          % Disturbance amplitude
omega_d = 1.5;    % Disturbance frequency

kv = 150;         % Velocity feedback gain

gamma = 1;        % Learning rate for RBF network

%% =================== Simulation Time ===================
dt = 0.01;
T = 10;
t = 0:dt:T;
N = length(t);

%% =================== Variable Initialization ===================
x = zeros(1,N);         % Position
v = zeros(1,N);         % Velocity
a = zeros(1,N);         % Acceleration
xd = zeros(1,N);        % Desired position
xd_dot = zeros(1,N);    % Desired velocity
xd_ddot = zeros(1,N);   % Desired acceleration
Fm = zeros(1,N);        % Control force
Fex = 5*sin(2*t);       % Human-robot interaction force
F_fric = zeros(1,N);    % Friction force
dist = zeros(1,N);      % Actual disturbance
d_hat = zeros(1,N);     % Estimated disturbance
e = zeros(1,N);         % Velocity tracking error

%% =================== RBF Network Initialization ===================
n_rbf = 30;                            % Number of RBF neurons
W_hat = zeros(n_rbf,1);               % Estimated weights
c = linspace(-3, 3, n_rbf);           % Centers of Gaussian functions
b_rbf = 1.5;                          % Width of Gaussian functions

%% =================== Main Simulation Loop ===================
for i = 2:N
    %% Step 1: Generate desired trajectory from admittance model
    xd_ddot(i) = (Fex(i) - b_d * xd_dot(i-1)) / m_d;
    xd_dot(i) = xd_dot(i-1) + xd_ddot(i) * dt;
    xd(i) = xd(i-1) + xd_dot(i) * dt;

    %% Step 2: Compute velocity tracking error
    e(i) = xd_dot(i) - v(i-1);

    %% Step 3: RBF disturbance estimation
    phi_z = exp(-((e(i) - c).^2) / (2 * b_rbf^2))';
    d_rbf = W_hat' * phi_z;
    d_hat(i) = d_rbf;

    %% Step 4: Actual friction and disturbance
    F_fric(i) = mu * sign(v(i-1)) + b * v(i-1);
    d_real = Ad * sin(omega_d * t(i));
    dist(i) = d_real;

    %% Step 5: Control law
    Fm(i) = m * (xd_ddot(i) + kv * e(i)) + F_fric(i) + d_hat(i) - Fex(i);

    %% Step 6: Update system dynamics
    a(i) = (Fm(i) + Fex(i) - F_fric(i) - d_real) / m;
    v(i) = v(i-1) + a(i) * dt;
    x(i) = x(i-1) + v(i) * dt;

    %% Step 7: Update RBF weights (modified)
    er = -m * (a(i) - (xd_ddot(i) + kv * e(i)));
    W_hat_dot = gamma * phi_z * er;
    W_hat = W_hat + W_hat_dot * dt;
end

%% =================== Plotting ===================
fontsize = 12;

figure('Position',[100,100,800,600]);

% ------ Figure 1: Desired vs Actual Velocity ------
subplot(4,1,1);
plot(t, xd_dot, 'r--', 'LineWidth', 1.5); hold on;
plot(t, v, 'b-', 'LineWidth', 1.5);
h1 = legend({'$\dot{x}_d$ (Desired)', '$v$ (Actual)'}, 'Interpreter', 'latex', 'FontSize', fontsize, 'Location', 'best');
ylabel('Velocity [m/s]', 'FontSize', fontsize, 'Interpreter', 'latex');
xlabel('Time [s]', 'FontSize', fontsize, 'Interpreter', 'latex');
title('Comparison of Desired and Actual Velocity', 'FontSize', fontsize, 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize);
grid on;

% ------ Figure 2: Velocity Tracking Error ------
subplot(4,1,2);
plot(t, e, 'm', 'LineWidth', 1.5);
ylabel('Error [m/s]', 'FontSize', fontsize, 'Interpreter', 'latex');
xlabel('Time [s]', 'FontSize', fontsize, 'Interpreter', 'latex');
title('Velocity Tracking Error $e(t)$', 'FontSize', fontsize, 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize);
grid on;

% ------ Figure 3: Disturbance and Friction Force ------
subplot(4,1,3);
plot(t, dist, 'k', 'LineWidth', 1.5); hold on;
plot(t, F_fric, 'r', 'LineWidth', 1.5);
h2 = legend({'$d(t)$ (Disturbance)', '$F_{\mathrm{fric}}$ (Friction)'}, 'Interpreter', 'latex', 'FontSize', fontsize, 'Location', 'best');
ylabel('Force [N]', 'FontSize', fontsize, 'Interpreter', 'latex');
xlabel('Time [s]', 'FontSize', fontsize, 'Interpreter', 'latex');
title('Disturbance and Friction Force', 'FontSize', fontsize, 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize);
grid on;

% ------ Figure 4: Motor Output Force ------
subplot(4,1,4);
plot(t, Fm, 'b', 'LineWidth', 1.5);
xlabel('Time [s]', 'FontSize', fontsize, 'Interpreter', 'latex');
ylabel('$F_m$ [N]', 'FontSize', fontsize, 'Interpreter', 'latex');
title('Motor Output Force $F_m$', 'FontSize', fontsize, 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize);
ylim_max = max(abs(Fm));      % Automatically determine symmetric y-limits
ylim([-ylim_max, ylim_max]);   % Symmetric y-axis range
grid on;
