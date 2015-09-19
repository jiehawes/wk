close all;
clear all;
Colour = hsv;  % colours for plots

% Parameters for 2 element Windkessel Model
R = 0.95;  % mmHg*s/cm^3, systemic peripheral resistance
C = 1.47;  % cm^3/mmHg, systemic arterial compliance

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

% figure(1);
% % Analysis over 'cycle' number of cardiac cycles
% for n = 1:cycle
%   t = (n - 1) * Tc : 0.01 : n * Tc;
%   % Blood flow for each cardiac cycle
%   Q = @(t) Q_0 * sine(t - (n - 1) * Tc) .* (t <= ((n - 1) * Tc + Ts));
%   subplot (3, 1, 1);
%   plot(t, Q(t), 'm-', 'LineWidth' , 2);
%   hold on;
%   xlim ([0, n * Tc]);
%   ylim ([0, 600]);
%   title('Aortic Blood Flow Model');
%   ylabel('Blood Flow (ml/s)');
%   xlabel('time (s)');
%   legend('Blood flow rate', 'Location', 'northeastoutside');
%   % Initial conditions all models
%   if (n == 1)
%     P_ss = 0;
%     P_ss2 = 0;
%   end
%
%   %% Analytical solution
%   % Analytical solution for systolic cycle
%   ts = (n - 1) * Tc : 0.01 : (n - 1) * Tc + Ts;
%   P0 = P_ss + Q_0 * Ts * R * (C * pi * R) / ((Ts^2 + C^2 * pi^2 * R^2));
%   P_s = @(t) P0 * exp(-((t - ts(1)) / (R *C))) - Q_0 * Ts * R * (C * pi * R *...
%       cos(pi * (t - ts(1)) / Ts) - Ts * sin(pi * (t - ts(1)) / Ts)) / (Ts^2 +...
%       C^2 * pi^2 * R^2);
%   P_sd = P_s(ts(end)); % IC for the diastolic phase
%
%   % Analytical solution for diastolic cycle
%   td = (n - 1) * Tc + Ts : 0.01 : n * Tc;
%   P_d = @(t) P_sd * exp(-(t - td(1)) / (R * C));
%   P_ss = P_d(td(end));  % IC for the systolic phase
%   subplot(3, 1, 2);
%   hold on;
%   plot(ts, P_s(ts), 'r', 'LineWidth', 2);
%   plot(td, P_d(td), 'b', 'LineWidth', 2);
%   ylim([0, 150]);
%   xlim([0, n * Tc]);
%   title('Aortic blood pressure (analytical solution)');
%   ylabel('Pressure(mmHg)');
%   xlabel('time (s)');
%   legend('Systolic Pressure', 'Diastolic Pressure', 'Location', 'northeastoutside');
%
%   %% Numerical Solution for 2 Element WM
%   t = (n - 1) * Tc : 0.01 : n * Tc;
%   Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
%       (t <= ((n - 1) * Tc + Ts));
%   Y2 = @(t, y2) (-y2 / (R * C) + Q(t) / C);
%   [t_m2, P_m2] = ode113(Y2, [(n - 1) * Tc; n * Tc], P_ss2, options);
%   P_ss2 = P_m2(end);
%   subplot(3, 1, 3);
%   hold on;
%   plot(t_m2, P_m2, 'g-', 'LineWidth', 2);
%   ylim([0, 150]);
%   xlim([0, cycle * Tc]);
%   title('Aortic blood pressure (numerical solution)');
%   ylabel('Pressure (mmHg)');
%   xlabel('time (s)');
%   legend('Numerical solution', 'Location', 'northeastoutside');
%
%   %% Extracting the Blood pressure values for all model for one cycle
%   if (n == 1)
%     t2s = ts;
%     P2s = P_s(ts); % Analytical solution for Systole
%     t2d = td;
%     P2d = P_d(td); % Analytical solution for Diastole
%     t2 = t_m2;
%     P2 = P_m2; % Numerical solution for 2 element WM
%   end
% end

%% Analysis of Model with Various Initial Conditions
% Effect of varying C for 2 element WM
% figure(1);
% for i = 1:3
%   for n = 1:cycle
%     if (n == 1)
%       P_ss = 0;
%     end
%     multiple = 1.0 + 0.2 * (i - 2);
%     t_m = (n - 1) * Tc : 0.01 : n * Tc;
%     Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
%         (t <= ((n - 1) * Tc + Ts));
%     Y_C = @(t, y2) (-y2 / (1.0 * R * multiple * C) + Q(t) / C / multiple);
%     [t_m, P_m] = ode113(Y_C, [(n - 1) * Tc; n * Tc], P_ss, options);
%     P_ss = P_m(end);
%     if (i == 1)
%       if (n == 1)
%         P1 = P_m;
%         T1 = t_m;
%       else
%         P1 = [P1; P_m];
%         T1 = [T1; t_m];
%       end
%     elseif (i == 2)
%       if (n == 1)
%         P2 = P_m;
%         T2 = t_m;
%       else
%         P2 = [P2; P_m];
%         T2 = [T2; t_m];
%       end
%     else
%       if (n == 1)
%         P3 = P_m;
%         T3 = t_m;
%       else
%         P3 = [P3; P_m];
%         T3 = [T3; t_m];
%       end
%     end
%     % stuff(i) = [stuff(i), P_m];
%     % plot(t_m, P_m, 'Color', Colour(i * 5, :));
%     % hold on;
%     % ylim([0, 200]);
%     % xlim([0, cycle * Tc]);
%     % title('Aortic Blood Pressure with varying C - 2 Element Windkessel');
%     % ylabel('Pressure (mmHg)');
%     % xlabel('time (s)');
%   end
% end
% plot(T1, P1, 'b-', T2, P2, 'r-', T3, P3, 'g-');
% ylim([0, 200]);
% xlim([0, 17]);
% title('Aortic blood pressure at R^{0} with varying C');
% ylabel('Pressure (mmHg)');
% xlabel('Time (s)');
% legend('C^{-1}', 'C^0', 'C^{+1}');

% Effect of varying R for 2 element WM
% scaling = 1.2;
% figure(2);
% for i = 1:3
%   for n = 1:cycle
%     if (n == 1)
%       P_ss = 0;
%     end
%     multiple = 1.0 + 0.1 * (i - 1);
%     t_m = (n - 1) * Tc : 0.01 : n * Tc;
%     Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
%         (t <= ((n - 1) * Tc + Ts));
%     Y_R = @(t, y2) (-y2 / (R * multiple * C * scaling) + Q(t) / C / scaling);
%     [t_m, P_m] = ode113(Y_R, [(n - 1) * Tc; n * Tc], P_ss, options);
%     P_ss = P_m(end);
%     plot(t_m, P_m, 'Color', Colour(i * 5, :));
%     if (i == 1)
%       if (n == 1)
%         P1 = P_m;
%         T1 = t_m;
%       else
%         P1 = [P1; P_m];
%         T1 = [T1; t_m];
%       end
%     elseif (i == 2)
%       if (n == 1)
%         P2 = P_m;
%         T2 = t_m;
%       else
%         P2 = [P2; P_m];
%         T2 = [T2; t_m];
%       end
%     else
%       if (n == 1)
%         P3 = P_m;
%         T3 = t_m;
%       else
%         P3 = [P3; P_m];
%         T3 = [T3; t_m];
%       end
%     end
%     % hold on;
%     % ylim([0, 200]);
%     % xlim([0, cycle * Tc]);
%     % title('Aortic Blood Pressure with varying R - 2 Element Windkessel');
%     % ylabel('Pressure (mmHg)');
%     % xlabel('time (s)');
%   end
% end
% plot(T1, P1, 'b-', T2, P2, 'r-', T3, P3, 'g-');
% ylim([0, 200]);
% xlim([0, 17]);
% title('Aortic blood pressure at C^{+1} with varying R');
% ylabel('Pressure (mmHg)');
% xlabel('Time (s)');
% legend('R^{-1}', 'R^0', 'R^{+1}');

scaling = 1.2;
figure(1);
for i = 1:2
  for n = 1:cycle
    if (n == 1)
      P_ss = 0;
    end
    t_m = (n - 1) * Tc : 0.01 : n * Tc;
    Q = @(t) (scale_fac * (n == 1) + 1.0) * Q_0 * sine(t - (n - 1) * Tc) .*...
        (t <= ((n - 1) * Tc + Ts));
    if (i == 1)
      Y_R = @(t, y2) (-y2 / (R * 0.8 * C * 1.2) + Q(t) / C / 1.2);
    else
      Y_R = @(t, y2) (-y2 / (R * 1.2 * C * 0.8) + Q(t) / C / 0.8);
    end
    [t_m, P_m] = ode113(Y_R, [(n - 1) * Tc; n * Tc], P_ss, options);
    P_ss = P_m(end);
    plot(t_m, P_m, 'Color', Colour(i * 5, :));
    if (i == 1)
      if (n == 1)
        P1 = P_m;
        T1 = t_m;
      else
        P1 = [P1; P_m];
        T1 = [T1; t_m];
      end
    else
      if (n == 1)
        P2 = P_m;
        T2 = t_m;
      else
        P2 = [P2; P_m];
        T2 = [T2; t_m];
      end
    end
  end
end
plot(T1, P1, 'b-', T2, P2, 'r-');
ylim([0, 200]);
xlim([0, 17]);
title('Aortic blood pressure in healthy and unhealthy individuals');
ylabel('Pressure (mmHg)');
xlabel('Time (s)');
legend('Healthy', 'Unhealthy');
