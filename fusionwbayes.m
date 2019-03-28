% main script to estimate intervals of mitdb records with added noise using
% different fusion algorithms as well as an optional Bayesian Process.
%
% Algorithms used:  Simple Weighted Average
%                   Median
%                   Best SQI
%                   Sigma-Trimmed Mean Filter
%                   Alpha-Trimmed Mean Filter
%                   Kalman Filter
%                   Best Bayes 
%
% Note: the prefix bayes_ refers to a matrix with the Bayesian Process
% applied
%
% Latest update: sigma trimmed filter no longer discardng intervals



% signals taken from mitdb to be used
signals = {'100','101','102','103','104',...
    '105','106','107','108',...
    '109','111','112','113',...
    '114','115','116','116',...
    '118','119','121','122',...
    '123','124','200','201',...
    '202','203','205','207',...
    '208','209','210','212',...
    '213','214','215','217',...
    '219','220','221','222',...
    '223','228','230','231',...
    '232','233','234'};

% strings to call files of relevant SNR 
dB = {'','e24-em','e18-em','e12-em','e06-em','e00-em','e_6-em','e_12-em','e_18-em','e_24-em','e_30-em',...
    'e24-ma','e18-ma','e12-ma','e06-ma','e00-ma','e_6-ma','e_12-ma','e_18-ma','e_24-ma','e_30-ma',...
    'e24-bw','e18-bw','e12-bw','e06-bw','e00-bw','e_6-bw','e_12-bw','e_18-bw','e_24-bw','e_30-bw'}; % strings to call files of relevant SNR 

f = waitbar(0,'beginning'); % waitbar for iterations

mae = cell(31,1); % cell containing data on signals
mean_mae = zeros(31,9);

bayes_mae = cell(31,1); % cell containing data on signals
bayes_mean_mae = zeros(31,8);

coverage = cell(31,1); % cell containing coverage data
mean_coverage = zeros(31,4);

% Data points to form bland altman plots for each individual
% algorithm
comp_swa = [];
comp_swa_b = [];
comp_bestSQI  = [];
comp_bestSQI_b = [];
comp_median = [];
comp_med_b = [];
comp_sigma = [];
comp_sig_b = [];
comp_alpha = [];
comp_alpha_b = [];
comp_kalman  = [];
comp_kalman_b = [];
comp_arg = [];
comp_swa_p = [];
comp_ind = [];
reference = [];

for snr = 1:31 % run over 31 different SNR's
    MAE = zeros(length(signals),9); % holds relevant information of errors (Mean Absolute Error)  
    
    BAYES_MAE = zeros(length(signals),8); % holds relevant information of errors (Mean Absolute Error of Bayesian Process)
    
    COVERAGE = zeros(length(signals),4); % holds information on coverage at SNR's
    
    for counting = 1:48 % 48 different files in database:
        tic
        signal_name = signals{counting};    % use specific file
        N0 = 1; % start pt of reading
        waitbar((snr-1)/31,f,strcat(signal_name,dB{snr})) % waitbar to show progress
        
        %%% First Lead: Lead II

        [signalorig]=rdsamp(strcat(signal_name,dB{snr}),1,[],N0);
        ann = rdann(signal_name,'atr',[],[],N0);
        c = length(ann);

        RR_control = diff(ann(1:c));  % our control RR-intervals for comparison
        RR_control(RR_control<0) = 0;
        annDet = cell(2,3); % annotation files
        
        % per sample SQI's
        fsqiDet = zeros(length(signalorig),2); 
        ksqiDet = zeros(length(signalorig),2);
        sqiDet = zeros(length(signalorig),2);
        
        size_arr = zeros(1,7); % array containing number fo detections for different detectors to match sizes 
        
        gap = 50; % window size for matching detections
        
        for i = 1:2 % fir 2 leads find detected beats
            [signalorig,Fs,tm]=rdsamp(strcat(signal_name,dB{snr}),i);
            
            %%% Detector 1: SQRS
            sqrs(strcat(signal_name,dB{snr}),[],[],i);
            annDet{i,1}=rdann(strcat(signal_name,dB{snr}),'qrs');
            fsqiDet(:,i) = fsqi_calc(signalorig,Fs);
            ksqiDet(:,i) = ksqi_calc(signalorig);
            sqiDet(:,i) = ksqiDet(:,i).*fsqiDet(:,i);
            size_arr(i*3-2) = length(annDet{i,1});

            %%% Detector 2 GQRS
            gqrs(strcat(signal_name,dB{snr}),[],[],i);
            annDet{i,2}=rdann(strcat(signal_name,dB{snr}),'qrs');
            size_arr(i*3-1) = length(annDet{i,2});
            
            %%% Detector 3 WQRS
            wqrs(strcat(signal_name,dB{snr}),[],[],i);
            annDet{i,3}=rdann(strcat(signal_name,dB{snr}),'wqrs');
            size_arr(i*3) = length(annDet{i,3});
        end
        
        size_arr(7) = c;
        d = max(size_arr); % Length of largest array
        
        %%% R-R Intervals
        
        ann_pre = zeros(d,7);
        for i = 1:2 % Resize so all are equal length
            ann_pre(:,3*i-2) = resize(annDet{i,1},d);  
            ann_pre(:,3*i-1) = resize(annDet{i,2},d);
            ann_pre(:,3*i) = resize(annDet{i,3},d);
        end
        
        % Do checks before convertinf beats to intervals
        check = diff(ann_pre);  % create RR-intervals
        ann_ind = ann_pre;      % RR for comparison only. No preprocessing step
        ann_pre((check<Fs/2.5)+1) = 0;  % discard any unreasonably small intervals 
        ann_pre(:,7) = resize(ann,d);    % annotated RR
        ann_pre(:,8:13) = ann_ind(:,1:6); % RR for comparison

        % Match and arrange the found beats so that they are matching
        [match,count] = arrange_ann(ann_pre, gap);

        % Timestamps of annotated beats
        ann_set = match(:,7);
        
        % Timestamp for individual detector for comparison
        ind_match = match(:,8:13);
        
        % Timestamps of qrs detectors groups in beats
        match = match(:,1:6);
        
        % Prepare intervals and sqi's
        interval = zeros(count-1,6);
        ind_interval = zeros(count-1,6);
        avg_sqi = zeros(count-1,6);

        % Need agreement between beats
        int = 0;
        
        for j = 1:count
            for i = 1:6
                if match(j,i) == 0 % if no value in column of row
                    int = int + 1; % increment case 
                end
                if int >= 4 % if less than 3 beats in row
                    match(j,:) = zeros(6,1); % set row to zero, decide not enough agreement 
                end
            end
            int = 0;
        end 
        
        % Create intervals 
        [interval,ind_interval,avg_sqi,empty_rows] = create_interval(match,ind_match,count,sqiDet,Fs);
        
        % Ground truth
        RR_set = zeros(count-1,1);
        
        % Create ground truth intervals grouped at indices corresponding to detected interval
        count_zero = 0;
        for i = 1:count-1
            if sum(ann_set(i+1:end)) == 0 
                break
            elseif empty_rows
                count_zero = count_zero - 1;
            elseif ann_set(i) > 0
                RR_set(i) = RR_control(i-count_zero);
                count_zero = count_zero - 1;
            end
            count_zero = count_zero + 1;
        end

        % using NaN's to easier ignore any 0's
        interval(interval == 0) = NaN;
   
        % apply Bayesian Process to intervals
        [p_of,p_val,best_bayes] = bayes(interval);

        % Remove unwanted intervals by Bayes
        bayes_interval = p_of.*interval;
        % Remove unwanted SQI's by Bayes and renormalise SQI
        bayes_avg_sqi = p_of.*avg_sqi;
        bayes_avg_sqi = bayes_avg_sqi./sum(bayes_avg_sqi,2);
                
        % Simple Weighted Average (SWA)
        avg_instRR = sum(interval.*avg_sqi,2,'omitnan')./sum(avg_sqi,2,'omitnan');  % take simple weighted avg instantaneous RR interval
        bayes_avg_instRR = sum(bayes_interval.*bayes_avg_sqi,2,'omitnan')./sum(bayes_avg_sqi,2,'omitnan');  % take simple weighted avg instantaneous RR interval
        
        avg_instRR(avg_instRR == 0) = NaN;
        bayes_avg_instRR(bayes_avg_instRR == 0) = NaN;
 
        % SWA with Bayes Porbabilities as SQI
        bayes_avg_p_instRR = sum(bayes_interval.*p_val.*bayes_avg_sqi,2,'omitnan')./sum(p_val.*bayes_avg_sqi,2,'omitnan');
        bayes_avg_p_instRR(bayes_avg_p_instRR == 0) = NaN;
        
        % Best SQI
        [bestSQI,bestInd] = max(avg_sqi,[],2);
        [bayes_bestSQI,bayes_bestInd] = max(bayes_avg_sqi,[],2);
        
        bestInd = [1:length(bestInd)]'+(bestInd-1).*(length(bestInd));
        bayes_bestInd = [1:length(bayes_bestInd)]'+(bayes_bestInd-1).*(length(bayes_bestInd));
        
        bestRR = interval(bestInd);
        bayes_bestRR = interval(bayes_bestInd);
        
        % Median
        medRR = median(interval,2,'omitnan');
        bayes_medRR = median(bayes_interval,2,'omitnan');
        
        % Alpha Trimmed Mean Filter
        alphaRR = alpha_trim(interval,avg_sqi);
        bayes_alphaRR = alpha_trim(bayes_interval,bayes_avg_sqi);  
        
        % Sigma Trimmed Mean Filter
        sigmaRR = zeros(length(interval),1);
        bayes_sigmaRR = zeros(length(interval),1);
        
        for i = 1:length(interval)
            
            % Vector containing data on intervals and SWA
            inquestion = [interval(i,:) avg_instRR(i)];
            bayes_inquestion = [bayes_interval(i,:) bayes_avg_instRR(i)];
            
            % only if there are more intervals at index 1
            if length(bayes_inquestion) > 2
                % Get mean and Standard Deviation
                meanRR = avg_instRR(i);
                bayes_meanRR = bayes_avg_instRR(i);

                stdRR = std(inquestion(inquestion>0),'omitnan');
                bayes_stdRR = std(bayes_inquestion(bayes_inquestion>0),'omitnan');

                % Remove intervals outside of half of Standard Deviation from SWA and recalculate mean 
                inquestion = inquestion(inquestion>(meanRR-stdRR/2));
                inquestion = mean(inquestion(inquestion<(meanRR+stdRR/2)),'omitnan');

                bayes_inquestion = bayes_inquestion(bayes_inquestion>=(bayes_meanRR-bayes_stdRR/2));
                bayes_inquestion = mean(bayes_inquestion(bayes_inquestion<=(bayes_meanRR+bayes_stdRR/2)),'omitnan');

                sigmaRR(i) = inquestion; 
                bayes_sigmaRR(i) = bayes_inquestion; 
            end
        end
        
        % Kalman Filter
        kalmanRR = dy_kalman_inter(interval,avg_instRR);
        bayes_kalmanRR = dy_kalman_inter(bayes_interval,bayes_avg_instRR);
                      
        % Best Bayes method
        best_bayes(best_bayes == 0) = NaN;
        
        % Individual Detector Error and Coverage
        tempcomp = zeros(6,1);
        coverage_individ = zeros(1,6);
        coverage_ind = zeros(1,6);
        
        % Coverage is ((length of intervals missed)-(length of signal)) 
        %                           /length of signal        
        for i = 1:6
            coverage_individ(i) = sum(RR_set((ind_interval(:,i) == 0))); 
            coverage_ind(i) = (length(signalorig)-coverage_individ(i))/length(signalorig);
            tempcomp(i) = sum(abs(ind_interval(:,i)-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);
        end
        
        % Find Coverage of Fusion Process, will be same for all algorithms
        unusable = find(isnan(avg_instRR));
        unusable = [unusable ; find(avg_instRR==0)];
        covering= sum(RR_set(unusable),'omitnan');
                        
        cover = (length(signalorig)-covering)/length(signalorig);
        
        % Genie Method will be the detector with lowest error with at least
        % 3/4 coverage of fusion or best coverage for each record
        cond = min([(3/4)*cover max(coverage_ind)]);
        
        [m_ind,tester] = min(tempcomp(coverage_ind >= cond));
        
        sensor_low = coverage_ind(tempcomp == m_ind);
        
        RR_set(RR_set == 0) = NaN;
        
        % Mean Absolute Errors
        MAE(counting,1) = mean(tempcomp(1:3));  % Mean error of lead 1
        MAE(counting,2) = mean(tempcomp(4:6));  % Meanerror of lead 2
        MAE(counting,3) = m_ind;  % Genie method
        MAE(counting,4) = sum(abs(avg_instRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);  % Err SWA  
        MAE(counting,5) = sum(abs(bestRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);      % Err Best SQI
        MAE(counting,6) = sum(abs(medRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);       % Err Median
        MAE(counting,7) = sum(abs(sigmaRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);     % Err Sigma
        MAE(counting,8) = sum(abs(alphaRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);     % Err Alpha
        MAE(counting,9) = sum(abs(kalmanRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);    % Err Kalman

        BAYES_MAE(counting,1) = sum(abs(bayes_avg_instRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);  % Err SWA  
        BAYES_MAE(counting,2) = sum(abs(bayes_bestRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);      % Err Best SQI
        BAYES_MAE(counting,3) = sum(abs(bayes_medRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);       % Err Median
        BAYES_MAE(counting,4) = sum(abs(bayes_sigmaRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);     % Err Sigma
        BAYES_MAE(counting,5) = sum(abs(bayes_alphaRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);     % Err Alpha
        BAYES_MAE(counting,6) = sum(abs(bayes_kalmanRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);    % Err Kalman
        BAYES_MAE(counting,7) = sum(abs(best_bayes-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);        % Err Best Bayes
        BAYES_MAE(counting,8) = sum(abs(bayes_avg_p_instRR-RR_set),'omitnan')/(length(RR_set(RR_set>0))*Fs);% Err SWA with Bayes as SQI
                     
        % Coverage: Mean Lead 1, Mean Lead 2, Best Individual Coverage,
        % Fusion
        COVERAGE(counting,:) = [mean(coverage_ind(1:3)) mean(coverage_ind(4:6)) sensor_low(1) cover];
               
        % Display Progress and Errors
        disp(signal_name)
        disp(MAE(counting,:))
        disp(BAYES_MAE(counting,:))
        
        % Create Array of all Data Points For creating Bland-Altman Plots
        % at 0dB Muscle Artifact
        if snr == 16
            comp_swa = [comp_swa; avg_instRR];  % SWA
            comp_swa = [comp_swa_b; bayes_avg_instRR];
            comp_bestSQI  = [comp_bestSQI; bestRR]; % Best SQI
            comp_bestSQI_b  = [comp_bestSQI_b; bayes_bestRR]; 
            comp_median = [comp_median; medRR];     % Median
            comp_med_b = [comp_med_b; bayes_medRR];
            comp_sigma = [comp_sigma; sigmaRR];     % Sigma
            comp_sig_b = [comp_sig_b; bayes_sigmaRR];
            comp_alpha = [comp_alpha; alphaRR];     % Alpha
            comp_alpha_b = [comp_alpha_b; bayes_alphaRR];
            comp_kalman  = [comp_kalman; kalmanRR]; % Kalman
            comp_kalman_b = [comp_kalman_b; bayes_kalmanRR];
            comp_arg = [comp_arg; best_bayes];              % Best Bayes
            comp_swa_p = [comp_swa_p; bayes_avg_p_instRR];  % SWA with Bayes SQI
            comp_ind = [comp_ind; ind_interval(:,tester)];  % "Genie"
            reference = [reference; RR_set];                % Reference Intervals
        end
    
        toc
    end
    
    % cells with all error information for current SNR 
    mae{snr} = MAE;
    disp(mean(MAE))
    mean_mae(snr,:) = mean(MAE);    
    bayes_mae{snr} = BAYES_MAE;
    bayes_mean_mae(snr,:) = mean(BAYES_MAE);
    coverage{snr} = COVERAGE;
    mean_coverage(snr,:) = mean(COVERAGE);
    
end
close(f) % Close Waitbar

% Plot Bland-Altman Plots
bland_altman(comp_swa,reference)
title('SWA ');
bland_altman(comp_swa_b,reference)
title('Bayes SWA ');
bland_altman(comp_bestSQI,reference)
title('Best SQI ');
bland_altman(comp_bestSQI_b,reference)
title('Bayes Best SQI ');
bland_altman(comp_median,reference)
title('Median ');
bland_altman(comp_med_b,reference)
title('Bayes Median ');
bland_altman(comp_sigma,reference)
title('Sigma ');
bland_altman(comp_sig_b,reference)
title('Bayes Sigma ');
bland_altman(comp_alpha,reference)
title('Bayes Alpha ');
bland_altman(comp_alpha,reference)
title('Alpha ');
bland_altman(comp_kalman_b,reference)
title('Bayes Kalman ');
bland_altman(comp_kalman,reference)
title('Kalman ');
bland_altman(comp_arg,reference)
title('Best Bayes ');
bland_altman(comp_swa_p,reference)
title('SQI with probabilities ');
bland_altman(comp_ind,reference)
title('Pr-Fusion Error');

% Make Other Plots of error to SNR
plotterer;

        
        