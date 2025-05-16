clear variables; close all; clc;

Npred = 2; % No. of predators
Nprey = 506; % No. of prey
Ncases = 50; % No. of realisations for each (chi, Omega_P) pair
chi = 17;
Lx = 2.5; Ly = 2.5; % Domain size
r_angle = [40,98,105,180]';

% Inputs for the DBSCAN algorithm
neighDist = 5.0 * (2.0*0.0136);
min_clustSize = round(0.01*Nprey);

for r = 1:numel(r_angle)
    figure(r)
    box on;
    filepath = sprintf('Fig2_data/chi_%d/rangle_%d', chi, r_angle(r));
    Nclusters = [];

    for case_no = 1:Ncases
        filename = sprintf('case_%d_data.txt', case_no);
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
        [idx correpts] = dbscan(D, neighDist, min_clustSize, 'Distance', 'precomputed');

        idx_filt = [x, y, idx];
        idx_filt = idx_filt(idx ~= -1,:);
        nclust = numel(unique(idx_filt(:,end)));
        for ns = 1:nclust
            Nclusters = [Nclusters; size(idx_filt(idx_filt(:,3)==ns,:),1)];
        end
    end
    Nclusters = Nclusters./Nprey; % Normalising the size of cluster with initial prey numbers
    histogram(Nclusters, 40, 'Normalization','probability')

    xlabel('$\varphi$','Interpreter','latex');
    ylabel('p($\varphi$)','Interpreter','latex')
    tle=sprintf('CSD $\\chi$ = %d, $\\Omega_{P}$ = %.3f $\\pi$',chi,(r_angle(r)/180));
    xlim([0 1.05]);ylim([0 1]);
    title(tle,'Interpreter','latex')
    
    outfile = sprintf('Fig2b_part%d.pdf', r);
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    print(gcf, outfile, '-dpdf', '-r0', '-bestfit');
    close all;
end