% Plots the changes in flock speed and local order parameter associated 
%with a split-join-split event
clear variables; close all; clc;
L = 2.5; dia_prey = 2*0.0136;
clust_speed = load('Fig3_cluster_speed.txt', '-ascii');
clust_speed = clust_speed.*sqrt(L/dia_prey);

loc_order = load('Fig3_locOP.txt', '-ascii');
loc_order(:,1) = loc_order(:,1).*sqrt(L/dia_prey);

xlabel('Time');
yyaxis left
plot(loc_order(:,1), loc_order(:,2), '-k')

ylabel('\Psi')
yyaxis right
hold on
plot(clust_speed(:,1), clust_speed(:,2), '-r')
plot(clust_speed(:,1), clust_speed(:,3), '--b')
ylabel('v_f^*');

saveas(gcf, 'Fig3.eps', 'epsc');

close all;

%% Plot insets
L = 1.0;
T = [14.85, 133.70, 230.23, 282.19, 341.61, 415.89, 504.96]';
npred = 2;

for t = 1:numel(T)
    filename = sprintf('data_%.2f.txt',T(t));
    file_dat = load(fullfile('Insets',filename), '-ascii');
    
    per_shift = file_dat(:,2)<0.5*L;
    file_dat(per_shift,2) = file_dat(per_shift,2)+L;
    file_dat(:,2) = file_dat(:,2)-0.5*L;

    pred_pos = file_dat(1:npred,:);
    file_dat(1:npred,:)=[];
    
    part_idx = file_dat(file_dat(:,end)~=-1, :);
    
    Nclust = numel(unique(part_idx(:,end)));
    
    box on; hold on;
    
    if Nclust > 1
        df_clust1 = file_dat(file_dat(:,end)==1,:);
        df_clust2 = file_dat(file_dat(:,end)==2,:);
        df_noClust = file_dat(file_dat(:,end)==-1,:);
        
        scatter(df_clust1(:,2), df_clust1(:,3), 20, 'r', 'filled')
        scatter(df_clust2(:,2), df_clust2(:,3), 20, 'b', 'filled')
        scatter(df_noClust(:,2), df_noClust(:,3), 20, 'm', 'filled')
        
        quiver(df_clust1(1:4:end,2), df_clust1(1:4:end,3), df_clust1(1:4:end,5), df_clust1(1:4:end,6), 0.5, 'k');
        quiver(df_clust2(1:4:end,2), df_clust2(1:4:end,3), df_clust2(1:4:end,5), df_clust2(1:4:end,6), 0.5, 'k');
        quiver(df_noClust(:,2), df_noClust(:,3), df_noClust(:,5), df_noClust(:,6), 0.1, 'k');
    else 
        df_clust1 = file_dat(file_dat(:,end)==1,:);
        df_noClust = file_dat(file_dat(:,end)==-1,:);
        
        scatter(df_clust1(:,2), df_clust1(:,3), 20, 'g', 'filled')
        scatter(df_noClust(:,2), df_noClust(:,3), 20, 'm', 'filled')
        
        quiver(df_clust1(1:4:end,2), df_clust1(1:4:end,3), df_clust1(1:4:end,5), df_clust1(1:4:end,6), 0.5, 'k');
        quiver(df_noClust(:,2), df_noClust(:,3), df_noClust(:,5), df_noClust(:,6), 0.1, 'k');
    end

    scatter(pred_pos(1,2), pred_pos(1,3), 160, 'k', 'filled')
    scatter(pred_pos(2,2), pred_pos(2,3), 160, 'r', 'filled')
    xlim([0 1])
    ylim([0 1])
    axis square
    
    outfile = sprintf('Fig3_inset_%.2f.eps', T(t));

    saveas(gcf, outfile, 'epsc');
    close all;
end
