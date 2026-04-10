%% voltage_comparison.m
% Compare two voltage traces in variables:
%   nv_voltage_trace  [Nx2 double]
%   eu_voltage_trace  [Mx2 double]
%
% Each has columns [time, voltage].

% Storage for ISI's
all_ISI = [];
labels = [];

% Path to neuron recordings
recording_path = "neuron_recordings/";

% Names for each ion channel config
% configs = [
%     "na_k_leak"
%     "na_k_leak_ca"
%     "na_k_leak_lowca"
%     "na_k_leak_slowk"
%     "na_k_leak_ca_lowca"
%     "na_k_leak_ca_slowk"
%     "na_k_leak_lowca_slowk"
%     "na_k_leak_ca_lowca_slowk"
% ];

% Choose prefixes for files
prefix1 = "yale_";
% prefix2 = "euler_";
prefix2 = "nv_";
% cfg = "na_k_leak";

for c = 1:numel(configs)
    cfg = configs(c); 

    % Build the two filepaths
    file1 = recording_path + prefix1 + cfg + ".csv";
    file2 = recording_path + prefix2 + cfg + ".csv";

    % Print what you're comparing (so logs are self-documenting)
    disp("Comparing config: " + cfg);
    disp("  data1 <- " + file1);
    disp("  data2 <- " + file2);

    % Load
    if ~isfile(file1) || ~isfile(file2)
        disp("  SKIP (missing one or both files)");
        continue;
    end

    % file1 = recording_path + prefix2 + cfg + ".csv";
    % file2 = recording_path + prefix1 + cfg + ".csv";
% data1 = readmatrix("C:\Users\playtoe\Programs\neuron\neurovisor-current\Assets\CSV_Files\neuro_visor_recording_02-06-2026-11-43-24.csv");
% data1 = readmatrix("/Users/carlos/Desktop/neuron/numerical_stability_in_neurovisor/neuron_recordings/yale_na_k_leak_ca.csv");
% data2 = readmatrix("/Users/carlos/Desktop/neuron/numerical_stability_in_neurovisor/neuron_recordings/euler_na_k_leak_ca.csv");
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
    
    % RMSE / Strength of Signal
    % Diff between max of each spike
    
    edges = linspace(min(t2), max(t2), 5);
    fprintf('\nInterval RMSE\n');
    
    for q = 1:4
        if q < 4
            mask = t2 >= edges(q) & t2 < edges(q+1);
        else
            mask = t2 >= edges(q) & t2 <= edges(q+1);
        end
    
        qrmse = sqrt(mean(diffVec(mask).^2));
    
        fprintf('[%.1f–%.1f ms] RMSE = %.6f mV\n', ...
            edges(q), edges(q+1), qrmse);
    end
    
    % 4) Discover Peaks
    % Peak indices via derivative sign change
    thresh = 35;

    dv1 = diff(v1i);
    pk1 = find(dv1(1:end-1) >= 0 & dv1(2:end) < 0) + 1;
    
    dv2 = diff(v2);  
    pk2 = find(dv2(1:end-1) >= 0 & dv2(2:end) < 0) + 1;  % handles plateaus    

    % filter out tiny peaks less that 0
    
    pk1 = pk1(v1i(pk1) > thresh);
    pk2 = pk2(v2(pk2)  > thresh);
    
    sigTotal = max(v2) - min(v2);
    nrmseTotal = rmseVal / sigTotal;
    fprintf('RMSE = %.6f mV | nRMSE = %.6f (RMSE/max diff)\n', rmseVal, nrmseTotal);

    % Peak times + simple time-alignment comparison (pair by order)
    n = min(numel(pk1), numel(pk2));
    tpk1 = t2(pk1(1:n));
    tpk2 = t2(pk2(1:n));
    
    dt = tpk2 - tpk1;
    rmse_dt = sqrt(mean(dt.^2));
    max_abs_dt = max(abs(dt));
    
    % fprintf('\nPeak Timing Comparison\n');
    % fprintf('Peaks found: Yale=%d, Euler=%d, Compared=%d\n', numel(pk1), numel(pk2), n);
    % fprintf('RMSE(|Δt|) = %.6f ms\n', rmse_dt);
    % fprintf('Max |Δt|   = %.6f ms\n', max_abs_dt);
    
    %% Print first few peak pairs
    % kshow = min(10, n);
    % fprintf('\nFirst %d peak times (ms):\n', kshow);
    % for k = 1:kshow
    %     fprintf('#%d: t1=%.6f, t2=%.6f, Δt=%.6f\n', k, tpk1(k), tpk2(k), dt(k));
    % end

    % 5) Intraspike Interval

    % Compute ISIs for both traces
    isi1 = diff(tpk1);  % Yale ISIs (ms)
    isi2 = diff(tpk2);  % Euler ISIs (ms)

    % Store mean ISIs
    all_ISI = [all_ISI; mean(isi1), mean(isi2)];
    labels = [labels; cfg];

    % Comparse ISIs of each program
    isi_diff = isi2 - isi1;
    rmse_isi = sqrt(mean(isi_diff.^2));
    mean_abs_isi_diff = mean(abs(isi_diff));
    max_abs_isi_diff = max(abs(isi_diff));
    
    fprintf('\nISI Comparison:\n');
    fprintf('Mean ISI difference: %.3f ms\n', mean(isi_diff));
    fprintf('RMSE(ISI): %.3f ms\n', rmse_isi);
    fprintf('Max ΔISI: %.3f ms\n', max_abs_isi_diff);

    % Instantaneous firing rate (1000 / ISI for Hz)
    fr1 = 1000 ./ isi1;  % Hz
    fr2 = 1000 ./ isi2;  % Hz

    fprintf('\nMean Firing Frequencies:\n');
    fprintf('  Yale:  %.2f Hz\n', mean(fr1));
    fprintf('  Euler: %.2f Hz\n', mean(fr2));
    fprintf('  Difference: %.2f Hz\n', mean(fr2) - mean(fr1));
    
    %% 6) Plot
    %figure('Name',"Voltage Comparison: "+cfg,'Color','white');
    %plot(t2, v2, 'b-', 'LineWidth',1.2); hold on;
    %plot(t2, v1i,'r--','LineWidth',1.0);
    %xlabel('Time (ms)'); ylabel('Voltage (mV)');
    %title(cfg);
    %legend({'Yale','Euler'}, 'Location','best');
    %grid on;

    %figure('Name',"Voltage Comparison: "+cfg,'Color','white','Position',[100 100 1200 800]);
    
    % Subplot 1: Voltage traces

    %plot(t2, v2, 'b-', 'LineWidth',1.2); hold on;
    %plot(t2, v1i,'r--','LineWidth',1.0);
    %plot(tpk2, v2(pk2(1:n)), 'bo', 'MarkerSize',6, 'MarkerFaceColor','b');  % NEW: peak markers
    %plot(tpk1, v1i(pk1(1:n)), 'ro', 'MarkerSize',6, 'MarkerFaceColor','r');

    % Create table (For ISI)
    T = table(labels, all_ISI(:,1), all_ISI(:,2), ...
        'VariableNames', {'Config','Yale_ISI','Euler_ISI'});
    
    disp(T);
    
    % Plot comparison
    figure;
    bar(all_ISI);
    set(gca, 'XTickLabel', labels);
    xlabel('Ion Channel Config');
    ylabel('Mean ISI (ms)');
    legend({'Yale','Euler'});
    title('ISI Comparison Across Configurations');
    grid on;

    subplot(3,1,1);
    plot(t2, v2, 'b-', 'LineWidth',1.2); hold on;
    plot(t2, v1i,'r--','LineWidth',1.0);
    xlabel('Time (ms)'); ylabel('Voltage (mV)');
    title(cfg);
    legend({'Yale','Euler'}, 'Location','best');
    grid on;
    
    % Subplot 2: ISI comparison
    subplot(3,1,2);
    spike_intervals = 1:numel(isi1);
    plot(spike_intervals, isi1, 'ro-', 'LineWidth',1.5); hold on;
    plot(spike_intervals, isi2, 'bo-', 'LineWidth',1.5);
    ylabel('ISI (ms)');
    title('Intraspike Intervals');
        
    % Subplot 3: ISI difference over time
    subplot(3,1,3);
    plot(spike_intervals, isi_diff, 'ko-', 'LineWidth',1.5);
    yline(0, 'r--');
    ylabel('ΔISI (ms) [Euler - Yale]');
    title('ISI Drift');
    
    % fprintf('%.6f, %.6f', max(data2), min(data2));
   
    % Close current plot before computing the next plot and stats.
    uiwait(gcf);
    
end

% | euler_isi - yale_isi | / yale_isi
% Send Zach an organized list of parameters
% - as detailed as possible. 
% ion channel panel, 