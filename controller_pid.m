clc;
clear;

%% ======================= System Parameters =======================
m_d =1;       % Admittance mass
b_d = 6;      % Admittance damping

m =0.5;       % Actual system mass
mu = 0.8;     % Coulomb friction coefficient
b = 0.4;      % Viscous friction coefficient

A_d =10;      % Disturbance amplitude
w_d = 1.5;    % Disturbance frequency

%% =================== Initial PID Parameters ======================
Kp0 =80; Ki0 =120; Kd0 = 0.05;
alpha1 =40; alpha2 =60; alpha3=0.02;

%% ==================== Simulation Parameters ======================
T = 10;
dt = 0.01;
t = 0:dt:T;
N = length(t);

%% =================== Initialize Variables ========================
Fe = zeros(1,N);         % Interaction force
F_m = zeros(1,N);        % Motor control force
u = zeros(1,N);          % PID output

x_d_dot = zeros(1,N);    % Desired velocity
x_d_ddot = zeros(1,N);   % Desired acceleration
x_d = zeros(1,N);        % Desired position

x = zeros(1,N);          % Actual position
x_dot = zeros(1,N);      % Actual velocity
x_ddot = zeros(1,N);     % Actual acceleration

e = zeros(1,N);          % Velocity error
int_e = 0; last_e = 0;

%% ====================== Simulation Loop ========================
for i = 2:N
    % --- Step 1: External interaction force (applied by patient) ---
    Fe(i) =5* sin(2 * t(i));

    % --- Step 2: Admittance control generates desired velocity ---
    x_d_ddot(i) = (Fe(i) - b_d * x_d_dot(i-1)) / m_d;
    x_d_dot(i) = x_d_dot(i-1) + x_d_ddot(i) * dt;
    x_d(i) = x_d(i-1) + x_d_dot(i) * dt;

    % --- Step 3: PID Controller (target: velocity tracking) ---
    e(i) = x_d_dot(i) - x_dot(i-1);
    int_e = int_e + e(i) * dt;
    de = (e(i) - last_e) / dt;

    % Adaptive gain adjustment
    Kp = Kp0 + alpha1 * tanh(abs(e(i)));
    Ki = Ki0 + alpha2 * tanh(abs(int_e));
    Kd = Kd0 + alpha3 * tanh(abs(de));

    % PID control output
    u(i) = Kp * e(i) + Ki * int_e + Kd * de;

    % --- Step 4: Friction force ---
    friction =mu *sign(x_dot(i-1)) + b * x_dot(i-1);

    % --- Step 5: Disturbance (sine) ---
    d_t =A_d * sin(w_d * t(i));

    % --- Step 6: Control law computation Fm ---
    F_m(i) = m * x_d_ddot(i)- Fe(i)+friction+u(i) ;

    % --- Step 7: Actual system dynamics ---
    x_ddot(i) = (F_m(i) + Fe(i) -friction-d_t) / m;
    x_dot(i) = x_dot(i-1) + x_ddot(i) * dt;
    x(i) = x(i-1) + x_dot(i) * dt;

    last_e = e(i);
end

%% ======================== Plotting ===========================
fontsize = 12; % Unified font size

figure;
% Set the figure window size (width 800px, height 600px)
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]
subplot(4,1,1);
plot(t, x_d_dot, 'r--', 'LineWidth', 1.5); hold on;
plot(t, x_dot, 'b', 'LineWidth', 1.2);
legend({'$\dot{x}_d$ (Desired)', '$\dot{x}$ (Actual)'}, ...
    'Interpreter','latex','FontSize',fontsize,'Location','best');
ylabel('Velocity [m/s]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Comparison of Desired and Actual Velocity', 'Interpreter','latex','FontSize',fontsize);
set(gca, 'FontSize', fontsize);
grid on;

subplot(4,1,2);
plot(t, e, 'r', 'LineWidth', 1.2);
ylabel('Error [m/s]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Velocity Tracking Error $e(t)$', 'Interpreter','latex','FontSize',fontsize);
set(gca, 'FontSize', fontsize);
grid on;

subplot(4,1,3);
plot(t, -A_d * sin(w_d * t), 'm', 'LineWidth', 1.2); hold on;
plot(t, -mu * sign(x_dot) - b * x_dot, 'r', 'LineWidth', 1.2);
legend({'$d(t)$ (Disturbance)', '$F_{\mathrm{fric}}$ (Friction)'}, ...
    'Interpreter','latex','FontSize',fontsize,'Location','best');
ylabel('Force [N]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Disturbance and Friction Force', 'Interpreter','tex','FontSize',fontsize);
set(gca, 'FontSize', fontsize);
grid on;

subplot(4,1,4);
plot(t, F_m, 'k', 'LineWidth', 1.2);
ylabel('$F_m$ [N]', 'Interpreter','latex','FontSize',fontsize);
xlabel('Time [s]', 'Interpreter','latex','FontSize',fontsize);
title('Motor Output Force $F_m$', 'Interpreter','latex','FontSize',fontsize);
set(gca, 'FontSize', fontsize);

ylim_max = max(abs(F_m));      % Automatically get the max absolute value
ylim([-ylim_max, ylim_max]);   % Symmetric y-axis setting

grid on;
