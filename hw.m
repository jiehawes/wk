close all;
clear all;
Colour = hsv;  % colours for plots

% Parameters for 2 element Windkessel Model
R = 0.95;  % mmHg/cm^3/s, systemic peripheral resistance
C = 1.06;  % cm^3/mmHg, systemic arterial compliance

% Simulation options, refine step size for ODE solver to produce smoother graphs
options = odeset('Refine', 16);
% scaling factor to increase blood flow for first cycle to get 0 to 120 mmHg
scale_fac = 1;
% number of cardiac cycles for which WM is analysed
cycle = 20;

%% Assumptions
Tc = 60 / 72;  % s, period of cardia cycle, 72 beats per second
Ts = (2 / 5) * Tc;  % s, period of systole

%% Modelling blood flow to the aorta
% Blood flow in one cardiac cycle is 90mL.
syms ti q
Q0 = solve(90 - int(q * (sin(pi * ti / Ts)), ti, 0, Ts), q);
Q_0 = eval(Q0);
sine = @(t) sin(pi * t / Ts);
% cardiac output for one cycle
Q = @(t) Q_0 * sine(t) .* (t <= Ts);

figure(1);
% Analysis over 'cycle' number of cardiac cycles
for n = 1:cycle
  t = (n - 1) * Tc : 0.01 : n * Tc;
  % Blood flow for each cardiac cycle
  Q = @(t) Q_0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));
  subplot (4, 1, 1);
  plot(t, Q(t), 'LineWidth' , 2);
  hold on;
  xlim ([0, n * Tc]);
  ylim ([0, 600]);
  title('Aortic Blood Flow Model');
  ylabel('Blood Flow (ml/s)');
  xlabel('time (s)');
  % Initial conditions all models
  if (n == 1)
    P_ss = 80;
    P_ss2 = 0;
  end

  %% Analytical solution
  % Analytical solution for systolic cycle
  ts = (n - 1) * Tc : 0.01 : (n - 1) * Tc + Ts;
  P0 = P_ss + Q_0 * Ts * R * (C * pi * R) / ((Ts^2 + C^2 * pi^2 * R^2));
  P_s = @(t) P0 * exp(-((t - ts(1)) / (R *C))) - Q_0 * Ts * R * (C * pi * R *...
      cos(pi * (t - ts(1)) / Ts) - Ts * sin(pi * (t - ts(1)) / Ts)) / (Ts^2 +...
      C^2 * pi^2 * R^2);
  P_sd = P_s(ts(end)); % IC for the diastolic phase

  % Analytical solution for diastolic cycle
  td = (n - 1) * Tc + Ts : 0.01 : n * Tc;
  P_d = @(t) P_sd * exp(-(t - td(1)) / (R * C));
  P_ss = P_d(td(end));  % IC for the systolic phase
  subplot(4, 1, 2);
  hold on;
  plot(ts, P_s(ts), 'r', 'LineWidth', 2);
  plot(td, P_d(td), 'b', 'LineWidth', 2);
  ylim([0, 150]);
  xlim([0, n * Tc]);
  title('Blood Pressure (Analytical Solution - 2 Element WM)');
  ylabel('Pressure(mmHg)');
  xlabel('time (s)');
  legend('Systolic Pressure', 'Diastolic Pressure');

  %% Numerical Solution for 2 Element WM
  t = (n - 1) * Tc : 0.01 : n * Tc;
  Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
      (t <= ((n - 1) * Tc + Ts));
  Y2 = @(t, y2) (-y2 / (R * C) + Q(t) / C);
  [t_m2, P_m2] = ode113(Y2, [(n - 1) * Tc; n * Tc], P_ss2, options);
  P_ss2 = P_m2(end);
  subplot(4, 1, 3);
  hold on;
  plot(t_m2, P_m2, 'LineWidth', 2);
  ylim([0, 150]);
  xlim([0, cycle * Tc]);
  title('Aortic Blood Pressure (Numerical Analysis - 2 Element WM)');
  ylabel('Pressure (mmHg)');
  xlabel('time (s)');
  legend('Numerical solution');

  %% Extracting the Blood pressure values for all model for one cycle
  if (n == 1)
    t2s = ts;
    P2s = P_s(ts); % Analytical solution for Systole
    t2d = td;
    P2d = P_d(td); % Analytical solution for Diastole
    t2 = t_m2;
    P2 = P_m2; % Numerical solution for 2 element WM
  end
end

%% Analysis of Model with Various Initial Conditions
% Effect of varying C for 2 element WM
figure(3);
for i = 1:3
  for n = 1:cycle
    if (n == 1)
      P_ss = 0;
    end
    multiple = 1.0 + 0.1 * (i - 1);
    t_m = (n - 1) * Tc : 0.01 : n * Tc;
    Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
        (t <= ((n - 1) * Tc + Ts));
    Y_C = @(t, y2) (-y2 / (R * multiple * C) + Q(t) / C / multiple);
    [t_m, P_m] = ode113(Y_C, [(n - 1) * Tc; n * Tc], P_ss, options);
    P_ss = P_m(end);
    plot(t_m, P_m, 'Color', Colour(i * 5, :));
    hold on;
    ylim([0, 200]);
    xlim([0, cycle * Tc]);
    title('Aortic Blood Pressure with varying C - 2 Element Windkessel');
    ylabel('Pressure (mmHg)');
    xlabel('time (s)');
  end
end

% Effect of varying R for 2 element WM
figure(4);
for i = 1:3
  for n = 1:cycle
    if (n == 1)
      P_ss = 0;
    end
    multiple = 1.0 + 0.1 * (i - 1);
    t_m = (n - 1) * Tc : 0.01 : n * Tc;
    Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
        (t <= ((n - 1) * Tc + Ts));
    Y_R = @(t, y2) (-y2 / (R * multiple * C) + Q(t) / C);
    [t_m, P_m] = ode113(Y_R, [(n - 1) * Tc; n * Tc], P_ss, options);
    P_ss = P_m(end);
    plot(t_m, P_m, 'Color', Colour(i * 5, :));
    hold on;
    ylim([0, 200]);
    xlim([0, cycle * Tc]);
    title('Aortic Blood Pressure with varying R - 2 Element Windkessel');
    ylabel('Pressure (mmHg)');
    xlabel('time (s)');
  end
end
