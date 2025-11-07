%% voltage_comparison.m
% Compare two voltage traces in variables:
%   nv_voltage_trace  [Nx2 double]
%   yn_voltage_trace  [Mx2 double]
%
% Each has columns [time, voltage].

disp("Starting voltage comparison script...");
data1 = readmatrix("C:\Users\playtoe\Programs\neuron\neurovisor-current\Assets\CSV_Files\neuro_visor_recording_11-03-2025-03-38-42.csv");
data2 = readmatrix("C:\Users\playtoe\Programs\neuron\yale\originalChannels_soma.csv");


% 1) Extract columns
t_nv = data1(:,1);
v_nv = data1(:,2);

t_ya = data2(:,1);
v_ya = data2(:,2);

% 2) Interpolate
v_nv_interp = interp1(t_nv, v_nv, t_ya, 'linear', 'extrap');

% 3) Differences
diffVec = v_ya - v_nv_interp;
mseVal  = mean(diffVec.^2);
rmseVal = sqrt(mseVal);
maxDiff = max(abs(diffVec));

% 5) Print stats
if mseVal == 0
    disp("mseVal is equal to 0.")
end
fprintf('\n=== Comparison Stats ===\n');
fprintf('MSE    = %.6f mV\n', mseVal);
fprintf('RMSE   = %.6f mV\n', rmseVal);
fprintf('MaxDiff= %.6f mV\n', maxDiff);

% 4) Plot
figure('Name','Euler vs NeuroVISOR','Color','white');
subplot(2,1,1);
plot(t_ya, v_ya, 'b-', 'LineWidth',1.5); hold on;
plot(t_ya, v_nv_interp, 'r--', 'LineWidth',1.2);
xlabel('Time (ms)'); ylabel('Voltage (mV)');
title('Euler vs NeuroVISOR');
legend({'Euler','NeuroVISOR'}, 'Location','best');
grid on;

subplot(2,1,2);
plot(t_ya, diffVec, 'k-', 'LineWidth',1.2);
xlabel('Time (ms)'); ylabel('Voltage Diff (mV)');
title('Difference: Euler - Neurovisor');
grid on;

% fprintf('%.6f, %.6f', max(data2), min(data2));

disp("Done with voltage comparison script.");
