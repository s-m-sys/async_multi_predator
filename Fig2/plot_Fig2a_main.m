clear variables; close all; clc;

Npred = 2; % No. of predators
Nprey = 506; % No. of prey
Ncases = 50; % No. of realisations for each (chi, Omega_P) pair
chi = [17,24,31]';
r_angle = [(10:10:80),(90:2:100),(110:10:180)]';
Lx = 2.5; Ly = Lx; % Domain size

% Inputs for DBSCAN algorithm
neighDist = 5.0*(2.0*0.0136);
min_clustSize = round(0.01*Nprey);

figure(1)
box on; hold on;
for i1 = 1:numel(chi)
    Nclust_Diag = zeros(numel(r_angle),2);
    for i2 = 1:numel(r_angle)
        Nclusters = zeros(Ncases,1); % Saves the no. of clusters for all cases for each (cv,r_angle) pair
        filepath = sprintf('Fig2_data/chi_%d/rangle_%d', chi(i1), r_angle(i2));
        
        for case_no = 1:Ncases
            filename = sprintf('case_%d_data.txt',case_no);
            dat = load(fullfile(filepath, filename), '-ascii');

            dat(1:Npred,:) = []; % Removing the predator data
            dat = dat(dat(:,1) > 0,:); % Removing the dead agents' data
            
            % Accounting for periodicity of domain
            x = dat(:,2);
            y = dat(:,3);
            D = zeros(size(x,1));
            for i = 1:numel(x)
                for j = i+1:numel(x)
                    dx = x(i)-x(j);
                    dy = y(i)-y(j);
                    
                    if dx > Lx/2
                        dx = dx-Lx;
                    elseif dx < -Lx/2
                        dx = dx+Lx;
                    end
                    if dy > Ly/2
                        dy = dy-Ly;
                    elseif dy < -Ly/2
                        dy = dy+Ly;
                    end
                    D(i,j) = sqrt(dx^2 + dy^2);
                    D(j,i) = D(i,j);
                end
            end
            % Cluster classification using DBSCAN algorithm
            [idx, correpts] = dbscan(D, neighDist, min_clustSize, 'Distance', 'precomputed');
            
            idx_filt = [x, y, idx];

            % figure(2)
            % gscatter(idx_filt(:,1), idx_filt(:,2), idx_filt(:,3));
            % xlim([0, Lx]); ylim([0 Ly]); axis square;grid minor;
            % outfile = sprintf('rAngle_%d_case_%d.eps',r_angle(i2),case_no);
            % filepath = sprintf('gscatter_plots/minN_%d/cv_%d',min_clustSize,cv);
            % saveas(gcf,fullfile(filepath,outfile),'epsc');

            idx_filt = idx_filt(idx ~= -1,:);
            idxx = unique(idx_filt(:,end));
            clustSize = zeros(numel(idxx),1);
            for i = 1:numel(idxx)
                clustSize(i,1) = size(idx_filt(idx_filt(:,end)==idxx(i),:),1);
            end
            Nclusters(case_no,1) = max(clustSize);
        end
        Nclust_Diag(i2,1) = mean(Nclusters./Nprey); % Normalising the size of cluster with
        Nclust_Diag(i2,2) = std(Nclusters./Nprey);  % initial prey numbers
    end
    plot(r_angle.*(pi/180),Nclust_Diag(:,1))
end

figure(1)
xlabel('\Omega_p'); ylabel('\langle l_c \rangle');
legend('\chi=17','\chi=24','\chi=31')
xticks(0:pi/4:pi); xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
yticks(0.5:0.1:1);
saveas(gcf, 'Fig2a_main.eps', 'epsc');

