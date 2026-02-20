%% voltage_comparison.m
% Compare two voltage traces in variables:
%   nv_voltage_trace  [Nx2 double]
%   yn_voltage_trace  [Mx2 double]
%
% Each has columns [time, voltage].

disp("Starting voltage comparison script...");
data1 = readmatrix("C:\Users\playtoe\Programs\neuron\neurovisor-current\Assets\CSV_Files\neuro_visor_recording_02-06-2026-11-43-24.csv");
% data1 = readmatrix("C:\Users\playtoe\Desktop\subject_01_eyesclosed_after.csv");
data2 = readmatrix("C:\Users\playtoe\Programs\neuron\numerical_stability_in_neurovisor\neuron_recordings\euler_trace.csv");


% 1) Extract columns
t_nv = data1(:,1);
v_nv = data1(:,2);

t_ya = data2(:,1);
v_ya = data2(:,2);

% 2) Interpolate
v_nv_interp = interp1(t_nv, v_nv, t_ya, 'linear', 'extrap');

% 3) Differences & global stats
diffVec = v_ya - v_nv_interp;
rmseVal = sqrt(mean(diffVec.^2));
fprintf('\n=== Comparison Stats ===\n');
fprintf('RMSE (total) = %.6f mV\n', rmseVal);

% Quadrant RMSE + spike counts (4 equal time windows)
thresh = 0;  % adjust as needed
t_min = min(t_ya); t_max = max(t_ya);
edges = linspace(t_min, t_max, 5);  % 4 intervals

fprintf('\n=== Per-Interval Stats ===\n');
for q = 1:4
    mask = t_ya >= edges(q) & t_ya < edges(q+1);
    if q == 4, mask = t_ya >= edges(q) & t_ya <= edges(q+1); end  % include endpoint
    qRMSE     = sqrt(mean(diffVec(mask).^2));
    spikes_eu = sum(diff(v_ya(mask) > thresh) == 1);
    spikes_nv = sum(diff(v_nv_interp(mask) > thresh) == 1);
    fprintf('[%.1f–%.1f ms] RMSE=%.4f mV | Spikes(Euler)=%d | Spikes(Yale)=%d\n', ...
        edges(q), edges(q+1), qRMSE, spikes_eu, spikes_nv);
end
fprintf('Total | RMSE=%.4f mV | Spikes(Euler)=%d | Spikes(Yale)=%d\n', ...
    rmseVal, sum(diff(v_ya > thresh)==1), sum(diff(v_nv_interp > thresh)==1));

% 5) Print stats
if rmseVal == 0
    rdisp("mseVal is equal to 0.")
end
fprintf('\n=== Comparison Stats ===\n');
% fprintf('MSE    = %.6f mV\n', mseVal);
fprintf('RMSE   = %.6f mV\n', rmseVal);
% fprintf('MaxDiff= %.6f mV\n', maxDiff);

% 4) Plot
figure('Name','Voltage Comparison','Color','white');
subplot(2,1,1,'Position',[0.05 0.55 0.9 0.4]);
plot(t_ya, v_ya, 'b-', 'LineWidth',1.5); hold on;
plot(t_ya, v_nv_interp, 'r--', 'LineWidth',1.2);
xlabel('Time (ms)'); ylabel('Voltage (mV)');
title('Voltage Comparison');
legend({'Euler','Yale'}, 'Location','best');
grid on;

% fprintf('%.6f, %.6f', max(data2), min(data2));

disp("Done with voltage comparison script.");
