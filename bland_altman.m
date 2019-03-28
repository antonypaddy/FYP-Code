function [] = bland_altman(in,ref)
% Function to plot a bland allman plot of the error between input an a
% reference
% Input in = input intervals
% Input ref = reference intervals

% calculate percentiles
[upper,lower,mean_err] = perc(in,ref,360,95);
count = 0;

% only include valid points and count
for i = 1:length(in)
    if in(i) > 0 && ref(i) > 0 
        count = count + 1;
        point(count) = (-in(i) + ref(i))/360;
        x_axis(count) = ref(i)/360;
    end
end

% randomly sample points
ind = randsample(count,round(count/30));

figure

% plot points
plot(x_axis(ind),point(ind),'+');
hold on
ylim([-0.75, 0.75]);
xlim([0, 2.2]);
xlabel('Reference RR interval (secs)')
ylabel('Error of method (secs)')
% create lines at percentiles
yline(upper);
yline(lower);
yline(mean_err);
text(2,upper+0.03,'E_9_5')
text(2,lower-0.025,'E_5')
text(2,mean_err+0.025,'Mean')
hold off

end

