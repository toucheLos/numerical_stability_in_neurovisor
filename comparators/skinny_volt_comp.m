%% skinny_volt_comp.m
% Compare two voltage traces in variables:
%   nv_voltage_trace  [Nx2 double]
%   eu_voltage_trace  [Mx2 double]
%
% Each has columns [time, voltage].

% Storage for ISI's
all_ISI = [];

file1 = "C:\Users\playtoe\Programs\neuron\neurovisor-current\Assets\CSV_Files\neuro_visor_recording_04-25-2026-04-39-55.csv";
file2 = ".\neuron_recordings\yale_na_k_leak";

data1 = readmatrix(file1);
data2 = readmatrix(file2);

fprintf('\nStarting voltage comparison script...\n\n');

% 1) Extract columns
t1 = data1(:,1); v1 = data1(:,2);
t2 = data2(:,1); v2 = data2(:,2);

% 2) Interpolate
v1i = interp1(t1, v1, t2, 'linear', 'extrap');

% 3) Differences & global stats
diffVec = v2 - v1i;
rmseVal = sqrt(mean(diffVec.^2));

% Print RMSE if 0
if rmseVal == 0
    disp("rmseVal is equal to 0.")
end

% Suggestions:
% RMSE / Strength of Signal
% Diff between max of each spike

% 4) Discover Peaks
% Peak indices via derivative sign change
thresh = 35;

dv1 = diff(v1i);
pk1 = find(dv1(1:end-1) >= 0 & dv1(2:end) < 0) + 1;

dv2 = diff(v2);  
pk2 = find(dv2(1:end-1) >= 0 & dv2(2:end) < 0) + 1;  % handles plateaus   

% Peak times + simple time-alignment comparison (pair by order)
n = min(numel(pk1), numel(pk2));
tpk1 = t2(pk1(1:n));
tpk2 = t2(pk2(1:n));

% 5) Intraspike Interval

% Compute ISIs for both traces
isi1 = diff(tpk1);  % Yale ISIs (ms)
isi2 = diff(tpk2);  % Euler ISIs (ms)

% Get the difference between ISI's
isi_diff = isi2 - isi1;

% Store mean ISIs
all_ISI = [all_ISI; mean(isi1), mean(isi2)];

fprintf('\nISI Comparison:\n');
fprintf('Mean ISI difference: %.3f ms\n', mean(isi_diff));

% 6) Plot
figure('Name',"Voltage Comparison", 'Color','white');
plot(t2, v2, 'b-', 'LineWidth',1.2); hold on;
plot(t2, v1i,'r--','LineWidth',1.0);
xlabel('Time (ms)'); ylabel('Voltage (mV)');
title("Na K Leak");
legend({'Yale','NeuroVISOR'}, 'Location','best');
grid on;