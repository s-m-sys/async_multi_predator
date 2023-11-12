clear variables; close all; clc;

T = (150:1:180)';
n_pred = 2;
h = 5*0.0136*2;
L = 2.5;
beta_e=15.0;

%% Plot Fig3a main

pred_traj = zeros(numel(T), 4);
clust_centroid = zeros(numel(T),4);
for t = 1:numel(T)
    file_name = sprintf('interior_%.2f.txt', T(t));
    file_dat = load(fullfile('Fig3_data/Fig3a_main', file_name), 'ascii');
    
    pred_traj(t, 1:2) = file_dat(1, 2:3); %Saving predator trajectories
    pred_traj(t, 3:4) = file_dat(2, 2:3);
    
    %% Remove negative ids and predator ids
    file_dat(1:n_pred,:) = [];
    file_dat = file_dat(file_dat(:,1) > 0, :);
    
    %% Building distance matrix
    N = size(file_dat,1); % No. of live prey
    D = zeros(N);
    for p1 = 1:N
    
        p2 = (1:N)';
    
        % Implementing the periodic boundary condition
        x_dist = file_dat(p1,2) - file_dat(p2,2);
        cond = x_dist > 0.5*L;
        x_dist(cond==1) = x_dist(cond==1) - L;
        cond = x_dist < -0.5*L;
        x_dist(cond==1) = x_dist(cond==1) + L;
    
        y_dist = file_dat(p1,3) - file_dat(p2,3);
        cond = y_dist > 0.5*L;
        y_dist(cond==1) = y_dist(cond==1) - L;
        cond = y_dist < -0.5*L;
        y_dist(cond==1) = y_dist(cond==1) + L;
    
        dist = sqrt(x_dist.^2 + y_dist.^2);
    
        D(p1,:) = dist;
    end

    %% Clustering based on precomputed distance
    [idx, correpts] = dbscan(D, h, 10, 'Distance', 'precomputed');
    
    file_dat = [file_dat, idx];
    file_dat = file_dat(file_dat(:,end) ~= -1, :);
    
    for n_clust = 1:2
        filter_dat = file_dat(file_dat(:,end) == n_clust, :);
        clust_centroid(t,(n_clust-1)*2+1:n_clust*2) = mean(filter_dat(:,2:3));
    end
end

clust_centroid = clust_centroid./L;
pred_traj = pred_traj./L;
singleClust = numel(clust_centroid(isnan(clust_centroid(:,3))));

%clust_centre = [clust_centroid(1:singleClust,1:2), repmat([0 0.4470 0.7410], singleClust,1)];
%clust_centre = [clust_centre; clust_centroid(singleClust+1:end,1:2), repmat([0.4660 0.6740 0.1880], size(clust_centroid,1)-singleClust,1)];
%clust_centre = [clust_centre; clust_centroid(singleClust+1:end,3:4), repmat([0.6350 0.0780 0.1840], size(clust_centroid,1)-singleClust,1)];
T = (T-T(1))./sqrt(L/beta_e);

close all;
box on; hold on;
scatter(clust_centroid(1:4:singleClust,1), clust_centroid(1:4:singleClust,2), 75, T(1:4:singleClust,1), 'filled');
scatter(clust_centroid(singleClust+1:2:end,1), clust_centroid(singleClust+1:2:end,2), 35, T(singleClust+1:2:end,1), 'filled');
scatter(clust_centroid(singleClust+1:2:end,3), clust_centroid(singleClust+1:2:end,4), 20, T(singleClust+1:2:end,1), 'filled');

scatter(pred_traj(1:5:end,1), pred_traj(1:5:end,2), 200, T(1:5:end,1), 'filled', 'MarkerEdgeColor', 'k');
scatter(pred_traj(1:5:end,3), pred_traj(1:5:end,4), 200, T(1:5:end,1), 'filled', 'MarkerEdgeColor', 'r');

colorbar
colormap cool

xlim([0 1]); ylim([0 1]);
axis square
%legend('P_1', 'P_2', 'F_0', 'F_1', 'F_2');
saveas(gcf, 'Fig3a_main.eps', 'epsc')