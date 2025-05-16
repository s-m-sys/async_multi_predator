clear variables; close all; clc;

Npred = 2; % No. of predators
Nprey = 506; % No. of prey
Ncases = 50; % No. of realisations for each (chi, Omega_P) pair
prey_dia = 0.0272;
chi = 31;
r_angle = [90,110,180]';
case_no = [7,34,46]';
Lx = 2.5; Ly = Lx; % Domain size

% Inputs for DBSCAN algorithm
neighDist = 5.0*(2.0*0.0136);
min_clustSize = round(0.01*Nprey);

Nclust_Diag = zeros(numel(r_angle),2);
for r = 1:numel(r_angle)
    Nclusters = zeros(Ncases,1); % Saves the no. of clusters for all cases for each (cv,r_angle) pair
    filepath = sprintf('Fig2_data/chi_%d/rangle_%d', chi, r_angle(r));

    filename = sprintf('case_%d_data.txt',case_no(r));
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

    idx_filt = [x./prey_dia, y./prey_dia, idx];

    figure(1)
    box on;
    gscatter(idx_filt(:,1), idx_filt(:,2), idx_filt(:,3),[],'o',2,'doleg','off',"filled");
    xlim([0, Lx./prey_dia]); ylim([0 Ly./prey_dia]);
    axis square;
    xlabel('x'); ylabel('y');
    tle = sprintf('$\\chi$ = %d, $\\Omega_{P}$ = %.3f $\\pi$', chi, r_angle(r)/180);
    title(tle, 'Interpreter','latex')
    outfile = sprintf('Fig2a_inset%d.eps', r);
    saveas(gcf,outfile,'epsc');
    close all;
end


