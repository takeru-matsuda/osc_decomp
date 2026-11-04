data = load(fullfile(fileparts(mfilename('fullpath')),'EEGExampleData.mat'));
y = data.eeg - mean(data.eeg);
fs = data.fs;
MAX_OSC = 6;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = ...
    osc_decomp(y,fs,MAX_OSC,'quasi-newton_grad1','whittle');
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_tau2 = osc_param(K,3*K+1);
fprintf('EEG: %s, %s, %.0f-%.0f s, fs = %g Hz.\n', ...
    data.metadata.record,data.metadata.channel,data.metadata.start_seconds, ...
    data.metadata.stop_seconds_exclusive,fs);
fprintf('AIC selects K = %d from 1:%d.\n',K,MAX_OSC);
fprintf('Frequencies (Hz):');
fprintf(' %.2f',osc_f);
fprintf('\n');
if any(osc_f == 0)
    fprintf('Zero-frequency components represent non-oscillatory background.\n');
end
if K == MAX_OSC
    fprintf('The selected K reaches MAX_OSC; consider increasing MAX_OSC.\n');
end

figure('Name','EEG data and model selection','Color','w');
subplot(2,1,1);
plot((1:numel(y))/fs,y,'k-');
xlim([1 numel(y)]/fs);
xlabel('Time in segment (s)');
ylabel('EEG (\muV)');
title('Eyes-closed EEG: S001R02, Oz, 10-20 s (mean removed)');
subplot(2,1,2);
plot(1:MAX_OSC,osc_AIC-minAIC,'ko-');
hold on;
plot(K,0,'ro','MarkerFaceColor','r');
xticks(1:MAX_OSC);
xlabel('Number of oscillators K');
ylabel('\DeltaAIC');
title(sprintf('AIC selects K = %d',K));

osc_plot(osc_mean,osc_cov,fs,K,0.95);
set(gcf,'Name','EEG oscillation components','Color','w');
for k = 1:K
    subplot(K,1,k);
    title(sprintf('Component %d: %.2f Hz',k,osc_f(k)));
    if osc_f(k) == 0
        title(sprintf('Component %d: 0 Hz (non-oscillatory background)',k));
    end
    ylabel('\muV');
end
xlabel('Time in segment (s)');
sgtitle('Smoothed components: pointwise 95% intervals');

phase_samples = 1:min(numel(y),round(3*fs));
eeg_rng_state = rng;
rng(1,'twister');
osc_phase_plot(osc_phase(:,phase_samples,:),osc_mean(:,phase_samples,:), ...
    osc_cov(:,:,phase_samples,:),fs,K,0.95,1000);
rng(eeg_rng_state);
set(gcf,'Name','EEG phase estimates','Color','w');
for k = 1:K
    subplot(K,1,k);
    title(sprintf('Component %d: %.2f Hz',k,osc_f(k)));
    ylabel('Phase (rad)');
    if osc_f(k) == 0
        cla;
        axis off;
        text(0.5,0.5,sprintf('Component %d: 0 Hz; phase not interpreted',k), ...
            'Units','normalized','HorizontalAlignment','center');
    end
end
xlabel('Time in segment (s)');
sgtitle('Phase: pointwise 95% intervals, first 3 s');

osc_spectrum_plot(y,fs,osc_a,osc_f,osc_sigma2,osc_tau2);
set(gcf,'Name','EEG spectrum','Color','w');
xlabel('Frequency (Hz)');
ylabel('log_{10} power density (\muV^2 / (rad/sample))');
title('EEG periodogram and fitted oscillator spectra');
