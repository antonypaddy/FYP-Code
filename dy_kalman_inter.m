function [out] = dy_kalman_inter(inter,control)
% Function that recreates a Kalman Filter for given intervals
% Input inter = matrix of intervals
% Input control = control intervals form different fusion algorithm to
% approximate errors
% Output out = Kalman Filtered result
%
% Prefix sig_ refers to an error
% delt refers to calculation error

control(isnan(control)) = 0;
inter(isnan(inter)) = 0;

est = [inter control]; % 

sig_delt = 0;   % error 
min_delt = 1e-6;% error must be greater that 0

% First error will be maximum interval of 
sig_delt(1) = max(abs(diff(control)))^2;

theta = zeros(length(inter)+1,1);
count = 0;
sig_i = ones(1,6);

for t = 2:length(inter)

    if sum(est(t-1,1:6)) > 0
        k = [theta(2:end) control];
        test = inter(t-1,:);
        d_len = length(test(test>0));
        V = ones((d_len),1);    
        R = diag(sig_i(test>0));
        sig_delt = max([sig_delt, min_delt]); 
        A = sig_delt*inv(sig_delt+R)*V;
        A = A./sum(A);  
        theta(t) = theta(t-count-1)*(1-(A'*V))+A'*test(test>0)';
        sig_i(test>0) = sig_i(test>0) + (((test(test>0)-control(t-1)).^2-sig_i(test>0))./t);
        sig_delt = sig_delt + (((control(t)-control(t-1))^2)-sig_delt)/1;
        count = 0;
    else
        count = count + 1;
    end
    
end 
out = theta(2:end); 
out(out == 0) = NaN;

end

