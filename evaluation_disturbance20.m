clc; clear; close all;

%% ========== Simulation Parameters ==========
T = 10;
dt = 0.01;
t = 0:dt:T;
N = length(t);

%% ========== System Parameters ==========
m = 0.5; mu = 0.8; b = 0.4;
A_d = 20; w_d = 1.5;
m_d = 1; b_d = 6;

%% Controller Switches
run_pid = true;
run_smc = true;
run_rbf = true;

smooth_sign = @(x) sign(x);

%% ========== PID Controller ==========
if run_pid
    Kp0 = 80; Ki0 = 120; Kd0 = 0.05;
    alpha1 = 40; alpha2 = 60; alpha3 = 0.02;
    pid_Fe = zeros(1,N); pid_u = zeros(1,N); pid_Fm = zeros(1,N);
    pid_x = zeros(1,N); pid_x_dot = zeros(1,N); pid_x_ddot = zeros(1,N);
    pid_xd = zeros(1,N); pid_xd_dot = zeros(1,N); pid_xd_ddot = zeros(1,N);
    pid_e = zeros(1,N); int_e = 0; last_e = 0;

    for i = 2:N
        pid_Fe(i) = 5 * sin(2 * t(i));
        pid_xd_ddot(i) = (pid_Fe(i) - b_d * pid_xd_dot(i-1)) / m_d;
        pid_xd_dot(i) = pid_xd_dot(i-1) + pid_xd_ddot(i) * dt;
        pid_xd(i) = pid_xd(i-1) + pid_xd_dot(i) * dt;

        pid_e(i) = pid_xd_dot(i) - pid_x_dot(i-1);
        int_e = int_e + pid_e(i) * dt;
        de = (pid_e(i) - last_e) / dt;

        Kp = Kp0 + alpha1 * tanh(abs(pid_e(i)));
        Ki = Ki0 + alpha2 * tanh(abs(int_e));
        Kd = Kd0 + alpha3 * tanh(abs(de));

        pid_u(i) = Kp * pid_e(i) + Ki * int_e + Kd * de;
        friction = mu * smooth_sign(pid_x_dot(i-1)) + b * pid_x_dot(i-1);
        d_t = A_d * sin(w_d * t(i));

        pid_Fm(i) = m * pid_xd_ddot(i) - pid_Fe(i) + friction + pid_u(i);
        pid_x_ddot(i) = (pid_Fm(i) + pid_Fe(i) - friction - d_t) / m;
        pid_x_dot(i) = pid_x_dot(i-1) + pid_x_ddot(i) * dt;
        pid_x(i) = pid_x(i-1) + pid_x_dot(i) * dt;

        last_e = pid_e(i);
    end
end

%% ========== SMC Controller ==========
if run_smc
    lambda = 20; rho = 0.8; alpha = 5; phi = 0.5; k1 = 30; k2 = 15;
    sat = @(x) tanh(x);
    smc_Fe = zeros(1,N); smc_Fm = zeros(1,N);
    smc_x = zeros(1,N); smc_x_dot = zeros(1,N); smc_x_ddot = zeros(1,N);
    smc_xd = zeros(1,N); smc_xd_dot = zeros(1,N); smc_xd_ddot = zeros(1,N);
    smc_e = zeros(1,N); smc_s = zeros(1,N);
    int_e = 0; int_sat_s = 0; smc_F_fric = zeros(1,N);

    for i = 2:N
        smc_Fe(i) = 5 * sin(2 * t(i));
        smc_xd_ddot(i) = (smc_Fe(i) - b_d * smc_xd_dot(i-1)) / m_d;
        smc_xd_dot(i) = smc_xd_dot(i-1) + smc_xd_ddot(i) * dt;
        smc_xd(i) = smc_xd(i-1) + smc_xd_dot(i) * dt;

        smc_e(i) = smc_xd_dot(i) - smc_x_dot(i-1);
        int_e = int_e + smc_e(i) * dt;
        smc_s(i) = smc_e(i) + lambda * int_e + rho * abs(smc_e(i))^alpha * sign(smc_e(i));
        sat_s = sat(smc_s(i) / phi);
        int_sat_s = int_sat_s + sat_s * dt;

        d_t = A_d * sin(w_d * t(i));
        smc_F_fric(i) = mu * smooth_sign(smc_x_dot(i-1)) + b * smc_x_dot(i-1);

        denominator = 1 + rho * abs(smc_e(i))^(alpha - 1) * sat_s;
        u_eq = smc_xd_ddot(i) + lambda * smc_e(i) + ...
               rho * alpha * abs(smc_e(i))^(alpha - 1) * smc_e(i) * sat_s + ...
               k1 * sqrt(abs(smc_s(i))) * sat_s + ...
               k2 * int_sat_s;

        smc_Fm(i) = m * u_eq / denominator - smc_Fe(i) + smc_F_fric(i);
        smc_x_ddot(i) = (smc_Fm(i) + smc_Fe(i) - smc_F_fric(i) - d_t) / m;
        smc_x_dot(i) = smc_x_dot(i-1) + smc_x_ddot(i) * dt;
        smc_x(i) = smc_x(i-1) + smc_x_dot(i) * dt;
    end
end

%% ========== RBF Controller ==========
if run_rbf
    kv = 150; gamma = 1;
    n_rbf = 30; b_rbf = 1.5;
    rbf_W = zeros(n_rbf,1); c = linspace(-3, 3, n_rbf);

    rbf_x = zeros(1,N); rbf_v = zeros(1,N); rbf_a = zeros(1,N);
    rbf_xd = zeros(1,N); rbf_xd_dot = zeros(1,N); rbf_xd_ddot = zeros(1,N);
    rbf_Fm = zeros(1,N); rbf_Fex = 5 * sin(2*t);
    rbf_F_fric = zeros(1,N); rbf_dist = zeros(1,N); rbf_d_hat = zeros(1,N);
    rbf_e = zeros(1,N);

    for i = 2:N
        rbf_xd_ddot(i) = (rbf_Fex(i) - b_d * rbf_xd_dot(i-1)) / m_d;
        rbf_xd_dot(i) = rbf_xd_dot(i-1) + rbf_xd_ddot(i) * dt;
        rbf_xd(i) = rbf_xd(i-1) + rbf_xd_dot(i) * dt;

        rbf_e(i) = rbf_xd_dot(i) - rbf_v(i-1);
        phi_z = exp(-((rbf_e(i) - c).^2) / (2 * b_rbf^2))';
        d_rbf = rbf_W' * phi_z;
        rbf_d_hat(i) = d_rbf;

        friction = mu * smooth_sign(rbf_v(i-1)) + b * rbf_v(i-1);
        rbf_F_fric(i) = friction;
        d_real = A_d * sin(w_d * t(i));
        rbf_dist(i) = d_real;

        rbf_Fm(i) = m * (rbf_xd_ddot(i) + kv * rbf_e(i)) + friction + d_rbf - rbf_Fex(i);
        rbf_a(i) = (rbf_Fm(i) + rbf_Fex(i) - rbf_F_fric(i) - d_real) / m;
        rbf_v(i) = rbf_v(i-1) + rbf_a(i) * dt;
        rbf_x(i) = rbf_x(i-1) + rbf_v(i) * dt;

        er = -m * (rbf_a(i) - (rbf_xd_ddot(i) + kv * rbf_e(i)));
        rbf_W_dot = gamma * phi_z * er;
        rbf_W = rbf_W + rbf_W_dot * dt;
    end
end

%% ========== Plotting ==========
fontsize = 12; linewidth = 1.2;
figure;
set(gcf, 'Position', [100, 100, 800, 600]);

% Plot desired and actual velocities
subplot(3,1,1); hold on; grid on;
if run_pid, plot(t, pid_xd_dot, 'k:', 'LineWidth', linewidth); plot(t, pid_x_dot, 'b', 'LineWidth', linewidth); end
if run_smc, plot(t, smc_x_dot, 'r--', 'LineWidth', linewidth); end
if run_rbf, plot(t, rbf_v, 'Color', [0.466 0.674 0.188], 'LineStyle', '-.', 'LineWidth', linewidth); end
legend({'$\dot{x}_d$ (Desired)','PID','SMC','RBF'},'Interpreter','latex','FontSize',fontsize,'Location','best');
title('Comparison of Desired and Actual Velocity','FontSize',fontsize);
ylabel('Velocity [m/s]','FontSize',fontsize);
xlabel('Time [s]','FontSize',fontsize);
set(gca,'FontSize',fontsize);

% Plot velocity tracking error
subplot(3,1,2); hold on; grid on;
if run_pid, plot(t, pid_e, 'b', 'LineWidth', linewidth); end
if run_smc, plot(t, smc_e, 'r--', 'LineWidth', linewidth); end
if run_rbf, plot(t, rbf_e, 'Color', [0.466 0.674 0.188], 'LineStyle', '-.', 'LineWidth', linewidth); end
legend({'PID','SMC','RBF'},'FontSize',fontsize);
title('Velocity Tracking Error','FontSize',fontsize);
ylabel('Error [m/s]','FontSize',fontsize);
xlabel('Time [s]','FontSize',fontsize);
set(gca,'FontSize',fontsize);

% Plot motor output force
subplot(3,1,3); hold on; grid on;
if run_pid, plot(t, pid_Fm, 'b', 'LineWidth', linewidth); end
if run_smc, plot(t, smc_Fm, 'r--', 'LineWidth', linewidth); end
if run_rbf, plot(t, rbf_Fm, 'Color', [0.466 0.674 0.188], 'LineStyle', '-.', 'LineWidth', linewidth); end
legend({'PID','SMC','RBF'},'FontSize',fontsize,'Location','best');
title('Motor Output Force $F_m$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('Force [N]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);
set(gca, 'FontSize', fontsize);
ylim_all = [pid_Fm smc_Fm rbf_Fm];
ylim_max = max(abs(ylim_all(:)));
ylim([-ylim_max, ylim_max]);

%% ========== Performance Evaluation ==========
threshold = 0.05; duration = 0.5;
calc_metrics = @(err, t, threshold, duration) struct( ...
    'rmse', sqrt(mean(err.^2)), ...
    'maxerr', max(abs(err)), ...
    'settle', calc_settling_time(err, t, threshold, duration) ...
);

if run_pid, metric_pid = calc_metrics(pid_e, t, threshold, duration); end
if run_smc, metric_smc = calc_metrics(smc_e, t, threshold, duration); end
if run_rbf, metric_rbf = calc_metrics(rbf_e, t, threshold, duration); end

%% ========== Display Performance Metrics ==========
if run_pid
    fprintf('\n== PID Controller Performance ==\n');
    fprintf('RMSE: %.5f m/s\n', metric_pid.rmse);
    fprintf('Max Absolute Error: %.5f m/s\n', metric_pid.maxerr);
    if isnan(metric_pid.settle)
        fprintf('Settling Time: N/A\n');
    else
        fprintf('Settling Time: %.5f s\n', metric_pid.settle);
    end
end

if run_smc
    fprintf('\n== Super-Twisting SMC Performance ==\n');
    fprintf('RMSE: %.5f m/s\n', metric_smc.rmse);
    fprintf('Max Absolute Error: %.5f m/s\n', metric_smc.maxerr);
    if isnan(metric_smc.settle)
        fprintf('Settling Time: N/A\n');
    else
        fprintf('Settling Time: %.5f s\n', metric_smc.settle);
    end
end

if run_rbf
    fprintf('\n== RBF Neural Network Controller Performance ==\n');
    fprintf('RMSE: %.5f m/s\n', metric_rbf.rmse);
    fprintf('Max Absolute Error: %.5f m/s\n', metric_rbf.maxerr);
    if isnan(metric_rbf.settle)
        fprintf('Settling Time: N/A\n');
    else
        fprintf('Settling Time: %.5f s\n', metric_rbf.settle);
    end
end

%% ========== Settling Time Calculation ==========
function settle_time = calc_settling_time(err, t, threshold, duration)
    dt = t(2) - t(1);
    window = round(duration / dt);
    for i = 1:(length(err) - window + 1)
        window_err = abs(err(i:i+window-1));
        if all(window_err < threshold)
            settle_time = t(i);
            return;
        end
    end
    settle_time = NaN;
end
