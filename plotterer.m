% Script to plot the MAE of certain algorithms vs SNR and Coverage

% Get Coverage
for i = 1:31
cov(i,:) = mean(coverage{i});
end

figure

% Plot Electrode Artifact MAE
h = plot(mean_mae(fliplr([1:11]),[3 8 9]));
set(h, {'color'}, {[1 0.1 0.1]; [0.1 1 0.1]; [0.5 0.5 0.5]});
set(h, {'Marker'}, {'*'; 'd'; 's'});
hold on

h = plot(bayes_mean_mae(fliplr([1:11]),[5 6 7 8]),'--');
ylim manual;
ylim([0,0.15])
set(h, {'color'}, {[0.1 1 0.1] ; [0.5 0.5 0.5]; [0.8 0.6 0.4]; [0.1 0.1 1]})
set(h, {'Marker'}, { 'd'; 's'; '+'; 'o'});
legend('Pre-Fusion','Alpha-Trimmed','Kalman','Alpha Trimmed','Kalman','Best Bayes','SWA w/ prob.')
 
title('Electrode Artifact')

ylabel('Mean Absolute Error (s)')

xlabel('SNR (dB)')

xticks(1:11)

xticklabels({'-30' '-24' '-18' '-12' '-6' '0' '6' '12' '18' '24' 'clean'})

% Electrode Artifact Coverage
figure

plot(cov(fliplr([1 12:21]),:))

title('Coverage Electrode Artifact')

ylabel('Coverage')

xlabel('SNR (dB)')

legend('Lead 1','Lead 2','Best det','Fusion')

xticks(1:11)

xticklabels({'30' '-24' '-18' '-12' '-6' '0' '6' '12' '18' '24' 'clean'})

% Muscle Artifact MAE
figure

h = plot(mean_mae(fliplr([1 12:21]),[3 8 9]));
set(h, {'color'}, {[1 0.1 0.1]; [0.1 1 0.1]; [0.5 0.5 0.5]});
set(h, {'Marker'}, {'*'; 'd'; 's'});
hold on

h = plot(bayes_mean_mae(fliplr([1 12:21]),[5 6 7 8]),'--');
ylim manual;
ylim([0,0.15])
set(h, {'color'}, { [0.1 1 0.1] ; [0.5 0.5 0.5]; [0.8 0.6 0.4]; [0.1 0.1 1]})
set(h, {'Marker'}, { 'd'; 's'; '+';'o'});
legend('Pre-Fusion','Alpha-Trimmed','Kalman','Alpha Trimmed','Kalman','Best Bayes','SWA w/ prob.')

title('Muscle Artifact')

ylabel('Mean Absolute Error (s)')

xlabel('SNR (dB)')

xticks(1:11)

xticklabels({'-30' '-24' '-18' '-12' '-6' '0' '6' '12' '18' '24' 'clean'})

% Muscle Artifact Coverage
figure

plot(cov(fliplr([1 12:21]),:))

title('Coverage Muscle Artifact')

ylabel('Coverage')

xlabel('SNR (dB)')

legend('Lead 1','Lead 2','Best det','Fusion')

xticks(1:11)

xticklabels({'30' '-24' '-18' '-12' '-6' '0' '6' '12' '18' '24' 'clean'})

% Baseline Wander MAE
figure

h = plot(mean_mae(fliplr([1 22:31]),[3 8 9]));
set(h, {'color'}, {[1 0.1 0.1];[0.1 1 0.1]; [0.5 0.5 0.5]});
set(h, {'Marker'}, {'*'; 'd'; 's'});
hold on

h = plot(bayes_mean_mae(fliplr([1 22:31]),[5 6 7 8]),'--');
ylim manual;
ylim([0,0.15])
set(h, {'color'}, {[0.1 1 0.1] ; [0.5 0.5 0.5];[0.8 0.6 0.4];  [0.1 0.1 1]})
set(h, {'Marker'}, {'d'; 's'; '+'; 'o'});
legend('Pre-Fusion','Alpha-Trimmed','Kalman','Alpha Trimmed','Kalman','Best Bayes','SWA w/ prob.')

title('Baseline Wander')

ylabel('Mean Absolute Error (s)')

xlabel('SNR (dB)')

xticks(1:11)

xticklabels({'-30' '-24' '-18' '-12' '-6' '0' '6' '12' '18' '24' 'clean'})

% Baseline Wander Coverage
figure

plot(cov(fliplr([1 22:31]),:))

title('Coverage Baseline Wander')

ylabel('Coverage')

xlabel('SNR (dB)')

legend('Lead 1','Lead 2','Best det','Fusion')

xticks(1:11)

xticklabels({'30' '-24' '-18' '-12' '-6' '0' '6' '12' '18' '24' 'clean'})
