clear variables; close all; clc;

dat = load('Fig6_data.txt', '-ascii');

errorbar(dat(:,1), dat(:,2), dat(:,3), 'ko')
xlabel('\chi');
ylabel('v^*');
%axis square

%% Fitting
f = fit(dat(:,1),dat(:,2),'poly2');
dat_fit = f.p1 .* (dat(:,1).^2) + f.p2 .* dat(:,1) + f.p3;
hold on;
plot(dat(:,1), dat_fit, 'b-')
legend('Numerical data', 'Scaling relation')

saveas(gcf, 'Fig6.eps', 'epsc')