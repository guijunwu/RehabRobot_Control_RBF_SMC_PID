clc; clear;

%% ========== Parameter Settings ==========
m_d = 1;        % Admittance mass
b_d = 6;        % Admittance damping
m = 0.5;        % Actual mass
mu = 0.8;       % Coulomb friction coefficient
b = 0.4;        % Viscous friction coefficient
A_d = 10;       % Disturbance amplitude
w_d = 1.5;      % Disturbance frequency

lambda = 20;    % Sliding surface parameter
rho = 0.8;      % Terminal sliding term coefficient
alpha = 5;      % Terminal sliding exponent
phi = 0.5;      % Boundary layer thickness for smoothing
k1 = 30;        % Gain for ST term
k2 = 15;        % Gain for integral term

sat = @(x) tanh(x);   % Smooth sat function

%% ========== Simulation Parameters ==========
T = 10;
dt = 0.01;
t = 0:dt:T;
N = length(t);

%% ========== Initialization ==========
x = zeros(1,N);         % Actual position
x_dot = zeros(1,N);     % Actual velocity
x_ddot = zeros(1,N);    % Actual acceleration

x_d = zeros(1,N);       % Desired position
x_d_dot = zeros(1,N);   % Desired velocity
x_d_ddot = zeros(1,N);  % Desired acceleration

e = zeros(1,N);         % Velocity error
de = 0;
s = zeros(1,N);         % Sliding surface
int_sat_s = 0;          % Integral of sat(s)

F_ex = zeros(1,N);      % Human-robot interaction force
F_m = zeros(1,N);       % Control input
F_fric = zeros(1,N);    % Friction force
d_t = zeros(1,N);       % External disturbance

%% ========== Main Simulation Loop ==========
int_e = 0;  % Initialize integral term

for i = 2:N
    % Step 1: External force + Admittance trajectory generation
    F_ex(i) = 5 * sin(2 * t(i));
    x_d_ddot(i) = (F_ex(i) - b_d * x_d_dot(i-1)) / m_d;
    x_d_dot(i) = x_d_dot(i-1) + x_d_ddot(i) * dt;
    x_d(i) = x_d(i-1) + x_d_dot(i) * dt;

    % Step 2: Error and sliding surface (Equation 16)
    e(i) = x_d_dot(i) - x_dot(i-1);
    int_e = int_e + e(i) * dt;
    s(i) = e(i) + lambda * int_e + rho * abs(e(i))^alpha * sign(e(i));

    % Step 3: Compute smoothed sat(s)
    sat_s = sat(s(i) / phi);
    int_sat_s = int_sat_s + sat_s * dt;

    % Step 4: Disturbance and friction
    d_t(i) = A_d * sin(w_d * t(i));
    F_fric(i) = mu * sign(x_dot(i-1)) + b * x_dot(i-1);

    % Step 5: Control law (Equation 21)
    denominator = 1 + rho * abs(e(i))^(alpha - 1) * sat_s;
    u_eq = x_d_ddot(i) + lambda * e(i) + ...
           rho * alpha * abs(e(i))^(alpha - 1) * e(i) * sat_s + ...
           k1 * sqrt(abs(s(i))) * sat_s + ...
           k2 * int_sat_s;

    F_m(i) = m * u_eq / denominator - F_ex(i) + F_fric(i);

    % Step 6: System dynamics update
    x_ddot(i) = (F_m(i) + F_ex(i) - F_fric(i) - d_t(i)) / m;
    x_dot(i) = x_dot(i-1) + x_ddot(i) * dt;
    x(i) = x(i-1) + x_dot(i) * dt;
end

%% ========== Plotting ==========
fontsize = 12; % Unified font size
figure;

% Set figure window size: width 800px, height 600px
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

subplot(4,1,1);
plot(t, x_d_dot, 'r--', 'LineWidth', 1.5); hold on;
plot(t, x_dot, 'b', 'LineWidth', 1.2);
h1 = legend({'$\dot{x}_d$ (Desired)', '$\dot{x}$ (Actual)'}, 'Interpreter','latex');
ylabel('Velocity [m/s]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Desired vs Actual Velocity', 'Interpreter','latex','FontSize',fontsize); 
set(gca, 'FontSize', fontsize);
set(h1,'FontSize',fontsize);
grid on;

subplot(4,1,2);
plot(t, e, 'r', 'LineWidth', 1.2);
ylabel('Error [m/s]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Tracking Error $e(t)$', 'Interpreter','latex','FontSize',fontsize); 
set(gca, 'FontSize', fontsize);
grid on;
ylim([-0.08 0.08]);
yticks([-0.08  0  0.08]);

subplot(4,1,3);
plot(t, -d_t, 'm', 'LineWidth', 1.2); hold on;
plot(t, -F_fric, 'g', 'LineWidth', 1.2);
h2 = legend({'$d(t)$ (Disturbance)', '$F_{\mathrm{fric}}(t)$ (Friction)'}, 'Interpreter','latex');
ylabel('Force [N]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Disturbance \& Friction Forces', 'Interpreter','latex','FontSize',fontsize); 
set(gca, 'FontSize', fontsize);
set(h2,'FontSize',fontsize);
grid on;

subplot(4,1,4);
plot(t, F_m, 'k', 'LineWidth', 1.2);
ylabel('$F_m$ [N]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Motor Output Force $F_m$', 'Interpreter','latex','FontSize',fontsize); 
set(gca, 'FontSize', fontsize);
ylim_max = max(abs(F_m));      % Automatically determine symmetric y-limits
ylim([-ylim_max, ylim_max]);   % Symmetric y-axis range
grid on;
