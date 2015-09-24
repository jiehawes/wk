clear all;
clf;
Colour = hsv;

% Simulation options, refine step size for ODE solver to produce smoother graphs
options = odeset('Refine', 8);

% Defining modeling parameters for Windkessel Model
% parameters for 3 element
R1= 0.05;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2 = 0.9278;  % mmHg*sec/cm^3, systemic peripheral resistance
C = 1.0666;  % cm^3/mmHg, systemic arterial compliance
%% Assumptions
Tc = 60 / 72;  % 72 beats per second
Ts = (2 / 5) * Tc;  % systole period
cycle = 1;  % number of cardiac cycles for which WM is analysed
% Modelling blood flow to the aorta
syms ti q
Qmax = solve(90 - int(q * (sin(pi * ti / Ts)), ti, 0, Ts), q);
Q0 = eval(Qmax);
sine = @(t) sin(pi * t / Ts);
Q = @(t) Q0 * sine(t) .* (t <= Ts); % for one cycle
figure(1);
% Analysis over 'cycle' number of cardiac cycles
% Initial conditions all models
P_ss = 80;
for n = 1:cycle
  t = (n - 1) * Tc : 0.01 : n * Tc;
  % Blood flow for each cardiac cycle
  Q = @(t) Q0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));

  %% Numerical Solution for 3 Element WM
  Y = @(t, y) (-y / (R2 * C) + Q(t) * (R2 + R1) / (R2 * C) + R1);
  [t_m, P_m] = ode113(Y, [(n - 1) * Tc; n * Tc], P_ss, options);
  P_ss = P_m(end);
  
  % Storing results in a single vector
  if (n == 1)
    P = P_m;
    T = t_m;
  else
    P = [P; P_m];
    T = [T; t_m];
  end
end

plot(T, P, 'LineWidth', 2);
ylim([0, 150]);
xlim([0, n * Tc]);
title('Aortic Blood Pressure (Numerical 3 Element WM)');
ylabel('Pressure (mmHg)');
xlabel('time (s)');
% Write to .dat file
dlmwrite('wk.dat', [T P], 'delimiter', ' ');
dlmwrite('wk1.dat', [t_m P_m], 'delimiter', ' ');
xlswrite('wk.xlsx', [T P]);
