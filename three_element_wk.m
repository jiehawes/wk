clear all;
clf;
Colour = hsv;

% Simulation options, refine step size for ODE solver to produce smoother graphs
options = odeset('Refine', 8);

% Defining modeling parameters for Windkessel Model
% parameters for 3 element
R1_3= 0.05;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_3 = 0.9278;  % mmHg*sec/cm^3, systemic peripheral resistance
C_3 = 1.47;  % cm^3/mmHg, systemic arterial compliance
% parameters for 4 element
R1_4 = 0.071;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_4 = 0.8495;  % mmHg*sec/cm^3, systemic peripheral resistance
C_4 = 1.47;  % cm^3/mmHg, systemic arterial compliance
L_4 = 5.2e-4;  % mmHs*s^2/cm^3, inertance of blood flow
%% Assumptions
Tc = 60 / 72;  % 72 beats per second
Ts = (2 / 5) * Tc;  % systole period
cycle = 1;  % number of cardiac cycles for which WM is analysed
% Modelling blood flow to the aorta
syms ti q
Qmax = solve(90 - int(q * (sin(pi * ti / Ts)), ti, 0, Ts), q);
Q0 = eval(Qmax);
sine = @(t) sin(pi * t / Ts);
cosine = @(t) cos(pi * t / Ts);
Q = @(t) Q0 * sine(t) .* (t <= Ts); % for one cycle
figure(1);
% Analysis over 'cycle' number of cardiac cycles
% Initial conditions all models
P_ss = 80;
P_ss2 = 80;
for n = 1:cycle
  t = (n - 1) * Tc : 0.01 : n * Tc;
  % Blood flow for each cardiac cycle
  Q = @(t) Q0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));
  dQdt = @(t) Q0 * pi / Ts * cosine(t - (n - 1) * Tc) .*...
      (t <= ((n - 1) * Tc + Ts));
  d2Qdt2 = @(t) Q0 * pi * pi / Ts / Ts * sine(t - (n - 1) * Tc) .*...
      (t <= ((n - 1) * Tc + Ts));
  %% Numerical Solution for 3 Element WM
  dydt = @(t, y) (-y / (R2_3 * C_3) + Q(t) * (R2_3 + R1_3) / (R2_3 * C_3) +...
      R1_3 * dQdt(t));
  [t_m, P_m] = ode113(dydt, [(n - 1) * Tc; n * Tc], P_ss, options);
  P_ss = P_m(end);

  %% Numerical solution for 4 Element WM
  dydt2 = @(t, y) (-y / (R2_4 * C_4) + Q(t) * (R2_4 + R1_4) / (R2_4 * C_4) +...
      (R1_4 + L_4 / C_4 / R2_4) * dQdt(t) + L_4 * d2Qdt2(t));
  [t_m2, P_m2] = ode113(dydt2, [(n - 1) * Tc; n * Tc], P_ss2, options);
  P_ss2 = P_m2(end);

  % Storing results in a single vector
  if (n == 1)
    P = P_m;
    T = t_m;
    P2 = P_m2;
    T2 = t_m2;
  else
    P = [P; P_m];
    T = [T; t_m];
    P2 = [P2; P_m2];
    T2 = [T2; t_m2];
  end
end

plot(T, P, 'LineWidth', 2);
hold on;
plot(T2, P2, 'LineWidth', 2);
disp(P_ss2);
ylim([0, 150]);
xlim([0, n * Tc]);
title('Aortic Blood Pressure (Numerical 3 Element WM)');
ylabel('Pressure (mmHg)');
xlabel('time (s)');
% % Write to .dat file
% dlmwrite('wk.dat', [T P], 'delimiter', ' ');
% dlmwrite('wk1.dat', [t_m P_m], 'delimiter', ' ');
% xlswrite('wk.xlsx', [T P]);
