clear all;
clf;
Colour = hsv;

% Simulation options, refine step size for ODE solver to produce smoother graphs
options = odeset('Refine', 8);

% Defining modeling parameters for Windkessel Model
% parameters for 3 element
R1_3= 0.05;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_3 = 1.1;  % mmHg*sec/cm^3, systemic peripheral resistance
C_3 = 1.37;  % cm^3/mmHg, systemic arterial compliance
% parameters for 4 element
R1_4a = 0.065;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_4a = 1.025;  % mmHg*sec/cm^3, systemic peripheral resistance
C_4a = 1.47;  % cm^3/mmHg, systemic arterial compliance
L_4a = 5.1e-4;  % mmHs*s^2/cm^3, inertance of blood flow
% parameters for 4 element (RL parallel)
R1_4b = 0.090;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_4b = 0.9;  % mmHg*sec/cm^3, systemic peripheral resistance
C_4b = 1.44;  % cm^3/mmHg, systemic arterial compliance
L_4b = 0.018;  % mmHs*s^2/cm^3, inertance of blood flow
%% Assumptions
Tc = 60 / 72;  % 72 beats per second
Ts = (2 / 5) * Tc;  % systole period
cycle = 1;  % number of cardiac cycles for which WM is analysed
% Modelling blood flow to the aorta
syms ti q
Qmax = solve(70 - int(q * (sin(pi * ti / Ts)), ti, 0, Ts), q);
Q0 = eval(Qmax);
sine = @(t) sin(pi * t / Ts);
cosine = @(t) cos(pi * t / Ts);
Q = @(t) Q0 * sine(t) .* (t <= Ts); % for one cycle
figure(1);
% Analysis over 'cycle' number of cardiac cycles
% Initial conditions all models
P_ss = 79.9;
P_ss2 = 80;
P_ss3 = 80;
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
  dydt2 = @(t, y) (-y / (R2_4a * C_4a) + Q(t) * (R2_4a + R1_4a) / (R2_4a *...
      C_4a) + (R1_4a + L_4a / C_4a / R2_4a) * dQdt(t) + L_4a * d2Qdt2(t));
  [t_m2, P_m2] = ode113(dydt2, [(n - 1) * Tc; n * Tc], P_ss2, options);
  P_ss2 = P_m2(end);

  %% Numerical solution for 4 element WM (LR parallel)
  d2ydt2 = @(t, y) [y(2); R1_4b * d2Qdt2(t) + (R1_4b / C_4b / R2_4b +...
      R1_4b / C_4b) * dQdt(t) + R1_4b / C_4b / L_4b * Q(t) - (1 / C_4b /...
      R2_4b + R1_4b / L_4b) * y(2) - R1_4b / C_4b / L_4b / R2_4b * y(1)];
  [t_m3, P_m3] = ode45(d2ydt2, [(n - 1) * Tc; n * Tc], [P_ss2; 0], options);
  P_ss3 = P_m3(end, 1);

  % Storing results in a single vector
  if (n == 1)
    P = P_m;
    T = t_m;
    P2 = P_m2;
    T2 = t_m2;
    P3 = P_m3(:, 1);
    T3 = t_m3;
  else
    P = [P; P_m];
    T = [T; t_m];
    P2 = [P2; P_m2];
    T2 = [T2; t_m2];
    P3 = [P3; P_m3(:, 1)];
    T3 = [T3; t_m3];
  end
end

plot(T, P, 'LineWidth', 2);
hold on;
plot(T2, P2, 'LineWidth', 2);
% hold on;
% plot(T3, P3, 'LineWidth', 2);
disp(P_ss);
disp(max(P));
ylim([0, 150]);
xlim([0, n * Tc]);
title('Aortic Blood Pressure (Numerical 3 Element WM)');
ylabel('Pressure (mmHg)');
xlabel('time (s)');
% Write to .dat file
dlmwrite('wk.dat', [T P], 'delimiter', ' ');
dlmwrite('wk1.dat', [t_m P_m], 'delimiter', ' ');
dlmwrite('wk2.dat', [t_m2 P_m2], 'delimiter', ' ');
xlswrite('wk.xlsx', [T P]);
