function [interval,ind_interval,sqi,empty_rows] = create_interval(match,ind_match,count,sqiDet,Fs)
% Function to create intervals from mathcing peaks and a second set of
% intervals relating to intervals we will be comparing with. We start with
% peak indices and we then switch interval indices. We will also remove any
% intervals < 0.333 secs or > 1.7 seconds
% 
% Input match = mathcing peaks to make intervals
% Input ind_match = matching peaks for comparison 
% Input count = length of files
% Input sqiDet = per sample SQI
% Input Fs =  Sampling Rate
% Output interval = intervals for fusion
% Output ind_interval = intervals for comparison
% Output sqi = sqi's of intervals
% Output empty_rows = index of empty rows 
 
    interval = zeros(count-1,6);
    ind_interval = zeros(count-1,6);
    sqi = zeros(count-1,6);
    empty_rows = zeros(count,1);
    % only calculate intervals of consecutive beats on lead
    for j = 1:count-1 % for length of files
        for i = 1:6 % for each detector 
            test = match(j+1:end,i); % vector of all next peaks from detector
            test = test(test>match(j,i)+Fs/3); % only create interval that is at least 0.333 seconds long
            match(j+find(match(j+1:end,i)<match(j,i)+Fs/3),i) = 0; % set unfeasible intervals to zero if less than 0.333 sec
            % if there is peak for comparison create interval between
            % current peak and next
            if ind_match(j,i) > 0 
                ind_test = ind_match(j+1:end,i);
                if find(ind_test)
                    ind_test = find(ind_test);
                    ind_interval(j,i) = ind_match(ind_test(1)+j,i) - ind_match(j,i);
                end
            end
            % if not any in current peak index
            if sum(match(j,:)) == 0 
                empty_rows(j) = 1; % indicate empty row
            % if no peak then set interval to zero and sqi at point to zero
            elseif match(j,i) == 0 || isnan(match(j,i))
                interval(j,i) = 0;
                sqi(j,i) = 0;
            % if curr and nxt peak != 0 
            elseif ~isnan(match(j,i)) &&  match(j,i) > 0 && ~isempty(test) 
                interval(j,i) = test(1) - match(j,i);  % calculate interval
                % if interval smaller than 1.7 seconds assign sqi 
                if interval(j,i) < Fs*1.7 
                    % SQI is mean of sqi over interval
                    sqi(j,i) = mean(sqiDet((match(j,i):test(1)),idivide(int64(i),int64(3), 'ceil'))); 
                end
            end
        end

    end
    sqi(interval>1.7*Fs) = 0;
    interval(interval>1.7*Fs) = 0;
    sqi = sqi./sum(sqi,2,'omitnan'); % nromalise SQI
    
    

end

