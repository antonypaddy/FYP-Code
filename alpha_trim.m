function [out] = alpha_trim(inter,sqi)
% Function to produce an Alpha Trimmed Weighted Mean of the intervals this
% involves removing certain extreme values either side of the median
% values
% Input inter = Matrix of intervals
% Input sqi = matrix of associated SQI's
% Output out = Alpha Trimmed-Mean 

for i = 1:length(inter)
    quest = sort(inter(i,:)); % Sort values in index row
    quest = quest(quest>0); % Remove 0's
    if length(quest) > 3 % Only apply if more than 3 intervals in this index are detected
        nz = quest([1 end]);
        down = find(inter(i,:) == nz(1)); % Remove largest and lowest values
        up = find(inter(i,:) == nz(2)); 
        sqi(i,down(1)) = 0; % Set relevant SQI and interval to 0
        inter(i,down(1)) = 0;
        sqi(i,up(1)) = 0;
        inter(i,up(1)) = 0;
    end
         
end

% Recalculate SWA with remaining intervals
out = sum(inter.*sqi,2,'omitnan')./sum(sqi,2,'omitnan');  % take simple weighted avg instantaneous RR interval

out(out == 0) = NaN;

end