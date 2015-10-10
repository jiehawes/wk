clear all;
clf;
Colour = hsv;

% Simulation options, refine step size for ODE solver to produce smoother graphs
options = odeset('Refine', 8);

% Defining modeling parameters for Windkessel Model
% parameters for 3 element
R1_3= 0.05;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_3 = 0.95;  % mmHg*sec/cm^3, systemic peripheral resistance
C_3 = 1.37;  % cm^3/mmHg, systemic arterial compliance
% parameters for 4 element
R1_4s = 0.065;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_4s = 1.025;  % mmHg*sec/cm^3, systemic peripheral resistance
C_4s = 1.47;  % cm^3/mmHg, systemic arterial compliance
L_4s = 5.1e-4;  % mmHs*s^2/cm^3, inertance of blood flow
% parameters for 4 element (RL parallel)
R1_4p = 0.090;  % mmHg*sec/cm^3, characteristic impedance of aorta
R2_4p = 0.9;  % mmHg*sec/cm^3, systemic peripheral resistance
C_4p = 1.44;  % cm^3/mmHg, systemic arterial compliance
L_4p = 0.018;  % mmHs*s^2/cm^3, inertance of blood flow
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
P_ss_wk3 = 79.9;
P_ss_wk4s = 80;
P_ss3 = 80;
P_ss_wk3_r1 = P_ss_wk3 * 0.8;
P_ss_wk3_r2 = P_ss_wk3 * 1.2;
P_ss_wk4s_r1 = P_ss_wk4s * 0.8;
P_ss_wk4s_r2 = P_ss_wk4s * 1.2;
% mu scale factor
r1_scale = 0.8;
r2_scale = 1.2;

for n = 1:cycle
  t = (n - 1) * Tc : 0.01 : n * Tc;
  % Blood flow for each cardiac cycle
  Q = @(t) Q0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));
  dQdt = @(t) Q0 * pi / Ts * cosine(t - (n - 1) * Tc) .*...
      (t <= ((n - 1) * Tc + Ts));
  d2Qdt2 = @(t) Q0 * pi * pi / Ts / Ts * sine(t - (n - 1) * Tc) .*...
      (t <= ((n - 1) * Tc + Ts));
  %% Numerical Solution for 3 Element WM
  dydt_wk3 = @(t, y) (-y / (R2_3 * C_3) + Q(t) * (R2_3 + R1_3) / (R2_3 *...
      C_3) + R1_3 * dQdt(t));
  [t_m_wk3, P_m_wk3] = ode113(dydt_wk3, [(n - 1) * Tc; n * Tc], P_ss_wk3,...
      options);
  P_ss_wk3 = P_m_wk3(end);

  %% Numerical solution for 4 Element WM
  dydt_wk4s = @(t, y) (-y / (R2_4s * C_4s) + Q(t) * (R2_4s + R1_4s) /...
      (R2_4s * C_4s) + (R1_4s + L_4s / C_4s / R2_4s) * dQdt(t) + L_4s *...
      d2Qdt2(t));
  [t_m_wk4s, P_m_wk4s] = ode113(dydt_wk4s, [(n - 1) * Tc; n * Tc],...
      P_ss_wk4s, options);
  P_ss_wk4s = P_m_wk4s(end);

  %% Numerical solution for 4 element WM (LR parallel)
  d2ydt2 = @(t, y) [y(2); R1_4p * d2Qdt2(t) + (R1_4p / C_4p / R2_4p +...
      R1_4p / C_4p) * dQdt(t) + R1_4p / C_4p / L_4p * Q(t) - (1 / C_4p /...
      R2_4p + R1_4p / L_4p) * y(2) - R1_4p / C_4p / L_4p / R2_4p * y(1)];
  [t_m_wk4p, P_m_wk4p] = ode45(d2ydt2, [(n - 1) * Tc; n * Tc],...
      [P_ss_wk4s; 0], options);
  P_ss3 = P_m_wk4p(end, 1);
  
  %% Variations to 0.8 mu
  % WK3
  R1_3_r1 = R1_3 * r1_scale;
  R2_3_r1 = R2_3 * r1_scale;
  dydt_wk3_r1 = @(t, y) (-y / (R2_3_r1 * C_3) + Q(t) * (R2_3_r1 + R1_3_r1) /...
      (R2_3_r1 * C_3) + R1_3_r1 * dQdt(t));
  [t_m_wk3_r1, P_m_wk3_r1] = ode113(dydt_wk3_r1, [(n - 1) * Tc; n * Tc],...
      P_ss_wk3_r1, options);
  P_ss_wk3_r1 = P_m_wk3_r1(end);
  % WK4
  R1_4s_r1 = R1_4s * r1_scale;
  R2_4s_r1 = R2_4s * r1_scale;
  dydt_wk4s_r1 = @(t, y) (-y / (R2_4s_r1 * C_4s) + Q(t) * (R2_4s_r1 +...
      R1_4s_r1) / (R2_4s_r1 * C_4s) + (R1_4s_r1 + L_4s / C_4s / R2_4s_r1) *...
      dQdt(t) + L_4s * d2Qdt2(t));
  [t_m_wk4s_r1, P_m_wk4s_r1] = ode113(dydt_wk4s_r1, [(n - 1) * Tc; n * Tc],...
      P_ss_wk4s_r1, options);
  P_ss_wk4s_r1 = P_m_wk4s_r1(end);
  
  %% Variations to 1.2 mu
  % WK3
  R1_3_r2 = R1_3 * r2_scale;
  R2_3_r2 = R2_3 * r2_scale;
  dydt_wk3_r2 = @(t, y) (-y / (R2_3_r2 * C_3) + Q(t) * (R2_3_r2 + R1_3_r2) /...
      (R2_3_r2 * C_3) + R1_3_r2 * dQdt(t));
  [t_m_wk3_r2, P_m_wk3_r2] = ode113(dydt_wk3_r2, [(n - 1) * Tc; n * Tc],...
      P_ss_wk3_r2, options);
  P_ss_wk3_r2 = P_m_wk3_r2(end);
  % WK4
  R1_4s_r2 = R1_4s * r2_scale;
  R2_4s_r2 = R2_4s * r2_scale;
  dydt_wk4s_r2 = @(t, y) (-y / (R2_4s_r2 * C_4s) + Q(t) * (R2_4s_r2 +...
      R1_4s_r2) / (R2_4s_r2 * C_4s) + (R1_4s_r2 + L_4s / C_4s / R2_4s_r2) *...
      dQdt(t) + L_4s * d2Qdt2(t));
  [t_m_wk4s_r2, P_m_wk4s_r2] = ode113(dydt_wk4s_r2, [(n - 1) * Tc; n * Tc],...
      P_ss_wk4s_r2, options);
  P_ss_wk4s_r2 = P_m_wk4s_r2(end);

  % Storing results in a single vector
  if (n == 1)
    P_wk3 = P_m_wk3;
    T_wk3 = t_m_wk3;
    P_wk4s = P_m_wk4s;
    T_wk4s = t_m_wk4s;
    P_wk4p = P_m_wk4p(:, 1);
    T_wk4p = t_m_wk4p;
    P_wk3_r1 = P_m_wk3_r1;
    T_wk3_r1 = t_m_wk3_r1;
    P_wk4s_r1 = P_m_wk4s_r1;
    T_wk4s_r1 = t_m_wk4s_r1;
    P_wk3_r2 = P_m_wk3_r2;
    T_wk3_r2 = t_m_wk3_r2;
    P_wk4s_r2 = P_m_wk4s_r2;
    T_wk4s_r2 = t_m_wk4s_r2;
  else
    P_wk3 = [P_wk3; P_m_wk3];
    T_wk3 = [T_wk3; t_m_wk3];
    P_wk4s = [P_wk4s; P_m_wk4s];
    T_wk4s = [T_wk4s; t_m_wk4s];
    P_wk4p = [P_wk4p; P_m_wk4p(:, 1)];
    T_wk4p = [T_wk4p; t_m_wk4p];
    P_wk3_r1 = [P_wk3_r1; P_m_wk3_r1];
    T_wk3_r1 = [T_wk3_r1; t_m_wk3_r1];
    P_wk4s_r1 = [P_wk4s_r1; P_m_wk4s_r1];
    T_wk4s_r1 = [T_wk4s_r1; t_m_wk4s_r1];
    P_wk3_r2 = [P_wk3_r2; P_m_wk3_r2];
    T_wk3_r2 = [T_wk3_r2; t_m_wk3_r2];
    P_wk4s_r2 = [P_wk4s_r2; P_m_wk4s_r2];
    T_wk4s_r2 = [T_wk4s_r2; t_m_wk4s_r2];
  end
end

plot(T_wk3, P_wk3, T_wk3_r1, P_wk3_r1, T_wk3_r2, P_wk3_r2, 'LineWidth', 2);
hold on;
plot(T_wk4s, P_wk4s, T_wk4s_r1, P_wk4s_r1, T_wk4s_r2, P_wk4s_r2, 'LineWidth', 2);

disp(P_ss_wk3);
disp(max(P_wk3));
ylim([0, 150]);
xlim([0, n * Tc]);
title('Blood pressure vs time');
ylabel('Pressure (mmHg)');
xlabel('time (s)');
legend('WK3 n', 'WK3 l', 'WK3 h', 'WK4 n', 'WK4 l', 'WK4 h', 'Location',...
    'northeastoutside');
% Write to .dat file
dlmwrite('wk.dat', [T_wk3 P_wk3], 'delimiter', ' ');
dlmwrite('wk3.dat', [t_m_wk3 P_m_wk3], 'delimiter', ' ');
dlmwrite('wk4s.dat', [t_m_wk4s P_m_wk4s], 'delimiter', ' ');
dlmwrite('wk3_r1.dat', [t_m_wk3_r1 P_m_wk3_r1], 'delimiter', ' ');
dlmwrite('wk4s_r1.dat', [t_m_wk4s_r1 P_m_wk4s_r1], 'delimiter', ' ');
dlmwrite('wk3_r2.dat', [t_m_wk3_r2 P_m_wk3_r2], 'delimiter', ' ');
dlmwrite('wk4s_r2.dat', [t_m_wk4s_r2 P_m_wk4s_r2], 'delimiter', ' ');
xlswrite('wk.xlsx', [T_wk3 P_wk3]);
