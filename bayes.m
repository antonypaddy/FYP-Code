function [p_out,p_val,arg] = bayes(inter)
% Function to apply Bayesian Process to intervals
% Input inter = intervals
% Output p_out = 1's and 0's of intervals to keep or discard
% Output p_val = normalised probabilities of intervals
% Output arg = most likely interval at each index

    threshold = 0.05; % threshold probability 
    inter = inter';
    inter(isnan(inter)) = 0;
    
    % define a new standard deviation that will be tighter than the normal
    % for this case using chi^2 distribution for 40 values given
    chi = 66.766;
    nz = inter(inter>0);
    f = movmean(nz,[40 0]);
    f(f<150) = 200;
    sigma = movstd(nz,[40 0]);
    sigma_l = sigma.*(sqrt(39/chi));
    sigma_l(1:40) = 60;
    % Calculate probability
    p = (1./(sqrt(2*pi).*sigma_l)).*(exp((-1/2).*(abs(f-nz)/chi).^2));
    p_out = zeros(size(inter));
    p_out(inter>0) = p;
    p_out(isinf(p_out)) = 0.95; 
    p_out = p_out';
    p_out = p_out./sum(p_out,2);    % Normalise the probabilites in indices
    p_val = p_out;                  
    [~,argmax] = max(p_out,[],2);   % Calculate index of most likely cases
    p_out(p_out<=threshold) = 0;    % Remove cases below threshold 
    p_out(p_out>threshold) = 1;     % set others to 1 for removing intervals
    index = [0:6:6*length(inter)-1]'+argmax;    % arrange indices of most likely to be in matrix index form  
    arg = inter(index); % find Best Bayes
    
end

