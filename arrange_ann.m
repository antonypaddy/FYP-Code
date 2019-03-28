function [match,count] = arrange_ann(ann, gap)
% Function that arranges a matrix on timestamps in windows of size
% determined by a parameter that first creates a single vector with all of
% the peaks and then matches and arranges them from then
%
% Input gap = window size
% Input ann = matrix of annotations
% Output match = matched peaks 
% Output count = length of peaks

    match = zeros(length(ann)*6,13);

    % Vector of all intervals arranged in order
    comp = zeros(2,length(ann(:)));
    temp = ann(:);

    % Arrange timestamps in order and keep timestamp indices 
    for i = 1:length(comp)
        [comp(1,i),comp(2,i)] = min(temp); %1st row beats, 2nd row detecor index
        temp(comp(2,i)) = NaN;
    end
    rem = comp(1,:) == 0;
    comp(:,rem) = [];

    % Zero padding for sake of comparing previous and next timestamps
    pos = [comp(1,1) comp(1,:) 0 0; 0 comp(2,:) 0 0];
    count = 1;

    % Group timestamps to likely corresponding beat
    % the following matches up relevant beats to each other according
    % to a windowsize decided earlier 
    for i = 2:length(pos)-2
        % arrange according to detector indices 
        place = idivide(int64(pos(2,i)), int64(length(ann)), 'ceil');
        if match(count,place) == 0                                          % if there is no beat in area
            if abs(pos(1,i)-pos(1,i+1)) > abs(pos(1,i)-pos(1,i-1)) + gap-1  % if next beat not in window        
                match(count,place) = pos(1,i);                              % place in window window
                count = count+1;                                            % increment window index
            elseif abs(pos(1,i)-pos(1,i-1)) > gap-1 % if current beat not in window
                count = count+1;                    % increment window index
                match(count,place) = pos(1,i);      % place beat
            else
                match(count,place) = pos(1,i);      % else place beat
            end
        else
            count = count + 1;                      % else incremnt window
            match(count,place) = pos(1,i);          % place beat
        end
    % start next beat is likely to correspond to new beat
    end

    % index will be 1 too large
    count = count - 1;
    match = match(1:count,:);
    match = match(any(match,2),:); % remove any signals that have no values to clean up and shorten loops
    count = length(match);
end

