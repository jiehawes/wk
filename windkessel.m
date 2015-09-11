close all;
clear all;
Colour = hsv;
% Defining modeling parameters f o r Windkessel Model
% parameters for 2 element
R = 0.95000;  % systemic peripheral resistance in (mmHg/cm^3/sec)
C = 1.0666;  % systemic arterial compliance in (cm^3/mmHg)
% parameters for 3 element
R1= 0.05;
R2 = 0.9000;  % systemic peripheral resistance in (mmHg/cm^3/sec)
C = 1.0666;
%% Assumptions
Tc = 60 / 72;  % 72 beats per second
Ts = (2 / 5) * Tc;  % systole period
cycle = 5;  % number of cardiac cycles for which WM is analysed
% Modelling blood flow to the aorta
syms ti q
I0 = solve(90 - int(q * (sin(pi * ti / Ts)), ti, 0, Ts), q);
I_0 = eval(I0);
sine = @(t) sin(pi * t / Ts);
I = @(t) I_0 * sine(t) .* (t <= Ts); % for one cycle
figure(1);
% Analysis over 'cycle' number of cardiac cycles
for n = 1:cycle
  t = (n - 1) * Tc : 0.01 : n * Tc;
  % Blood flow for each cardiac cycle
  I = @(t) I_0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));
  subplot (4, 1, 1);
  plot(t, I(t), 'LineWidth' , 2);
  hold on;
  xlim ([0, n * Tc]);
  ylim ([0, 600]);
  title('Aortic Blood Flow Model');
  ylabel('Blood Flow (ml/s)');
  xlabel('time (s)');
  % Initial conditions all models
  if (n == 1)
    P_ss = 80;
    P_ss2 = 80;
    P_ss3 = 80;
  end

  %% Analytical solution

  % Analytical solution for systolic cycle
  ts = (n - 1) * Tc : 0.01 : (n - 1) * Tc + Ts;
  P0 = P_ss + I_0 * Ts * R * (C * pi * R) / ((Ts^2 + C^2 * pi^2 * R^2));
  P_s = @(t) P0 * exp(-((t - ts(1)) / (R *C))) - I_0 * Ts * R * (C * pi * R *...
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
  I = @(t) I_0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));
  Y2 = @(t, y2) (-y2 / (R * C) + I(t) / C);
  [t_m2, P_m2] = ode45(Y2, [(n - 1) * Tc; n * Tc], P_ss2);
  P_ss2 = P_m2(end);
  subplot(4, 1, 3);
  hold on;
  plot(t_m2, P_m2, 'LineWidth', 2);
  ylim([0, 150]);
  xlim([0, cycle * Tc]);
  title('Aortic Blood Pressure (Numerical Analysis - 2 Element WM)');
  ylabel('Pressure (mmHg)');
  xlabel('time (s)');

  %% Numerical Solution for 3 Element WM
  Y3 = @(t, y3) (-y3 / (R2 * C) + I(t) * (R2 + R1) / (R2 * C) + R1);
  [t_m3, P_m3] = ode45(Y3, [(n - 1) * Tc; n * Tc], P_ss3);
  P_ss3 = P_m3(end);
  subplot(4, 1, 4);
  plot(t_m3, P_m3, 'LineWidth', 2);
  hold on;
  ylim([0, 150]);
  xlim([0, n * Tc]);
  title('Aortic Blood Pressure (Numerical 3 Element WM)');
  ylabel('Pressure (mmHg)');
  xlabel('time (s)');

  %% Extracting the Blood pressure values for all model for one cycle
  if (n == 1)
    t2s = ts;
    P2s = P_s(ts); % Analytical s o l u t ion f o r Sys t o le
    t2d = td;
    P2d = P_d(td); % Analytical s o l u t ion f o r Diastole
    t2 = t_m2;
    P2 = P_m2; % Numerical s o l u t i on f o r 2 element WM
    t3 = t_m3;
    P3 = P_m3; % Numerical s o l u t i on f o r 3 element WM
  end
end

%% Analysis of Model with Various Initial Conditions
% Modelling Blood flow into Aorta
squared = @(t) (square(pi * (t - ((Tc / 2) - Ts)) / Ts, 100 * (2 * Ts -...
    (Tc / 2)) / Tc) + 1) / 2;
% Estimated perodic Blood flow model
I_estimate = @(t) I_0 * squared(t) .* sine(t);
% Effect of varying IC for 2 element WM
figure(3);
Y = @(t, y) (-y / (R * C) + I_estimate(t) / C);
for i = 1:10
  [t_m, P_m] = ode45(Y, [0; cycle * Tc], 30 + i * 10);
  plot(t_m, P_m, 'Color', Colour(i * 5, :));
  hold on;
  ylim([0, 200]);
  xlim([0, cycle * Tc]);
  title('Aortic Blood Pressure with varying IC- 2 Element Windkessel');
  ylabel('Pressure (mmHg)');
  xlabel('time (s)');
end
% Effect of varying IC for 2 element WM
figure(5)
Y = @(t, y) (-y / (R2 * C) + I_estimate(t) * (R2 + R1) / (R2 * C) + R1);
for i = 1:10
  [t_m, P_m] = ode45(Y, [0; cycle * Tc], 30 + i * 10);
  plot(t_m , P_m , 'Color' , Colour(i * 5, :));
  hold on;
  ylim([0, 200]);
  xlim([0, cycle * Tc]);
  title('Aortic Blood Pressure with varying IC (3 element)');
  ylabel('Pressure (mmHg)');
  xlabel('time (s)');
end

%% Comparison of WM and analytical Solutions
figure(4);
hold on;
plot(t2s, P2s, 'r -.*', t2d, P2d, 'b-o', t2, P2, 'm-s', t3, P3, 'g-x',...
    'LineWidth', 2, 'MarkerSize', 5);
legend ('2 element WM-Analytical (Systolic)',...
    '2 element WM-Analytical (Diastolic)', '2 element WM', '3 element WM');
ylim([0, 150]);
xlim([0, Tc]);
