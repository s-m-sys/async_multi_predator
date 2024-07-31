clear variables; close all; clc;

chi = [163, 229, 294]';

%% Plot Fig4d main
colcode = ['r','g','m'];
mkrcode = ['o','^','x'];

box on; hold on;
for i = 1:numel(chi)
    filename = sprintf('chi_%d_nclusters.txt', chi(i));
    dat = load(fullfile('Fig4_data/Fig4d_data',filename), '-ascii');
    
    scatter(dat(:,1), dat(:,2), [], colcode(i), mkrcode(i))
end

ylim([0, 1.8])
xlim([0 200])
saveas(gcf, 'Fig4d_main.eps', 'epsc')
close all;

%% Plot Fig4d inset
dat = load(fullfile('Fig4_data/Fig4d_data', 'nclusters_relang_100.txt'), '-ascii');
bars_clust = zeros(2, numel(chi)+1); %First row for herding, second row for splitting
for i = 1:numel(chi)
    bars_clust(:,1) = [2;1];
    bars_clust(1,i+1) = numel(dat(dat(:,i)==2, i))/size(dat,1);
    bars_clust(2,i+1) = numel(dat(dat(:,i)==1, i))/size(dat,1);
end
barh(bars_clust(:,1), bars_clust(:,2:end)) %1 for herding, 2 for splitting
saveas(gcf, 'Fig4d_inset.eps', 'epsc')