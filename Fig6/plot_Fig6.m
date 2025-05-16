clear variables; close all; clc;

%% Plot Fig 6a
figure(1)
box on; hold on;
for chi = [17,31]
    infile = sprintf('Np_1_chi_%d_v2.txt', chi);
    prey_data_1pred = load(infile, '-ascii');
    infile = sprintf('Np_2_chi_%d_v2.txt', chi);
    prey_data_2pred = load(infile, '-ascii');

    plot(prey_data_1pred(:,1), prey_data_1pred(:,2))
    plot(prey_data_2pred(:,1), prey_data_2pred(:,2))
end

figure(1)
xlabel('\tau'); ylabel('N_l');
legend('N_p=1, \chi=17', 'N_p=2, \chi=17', 'N_p=1, \chi=31', 'N_p=2, \chi=31');
title('Fig6a');
saveas(gcf, 'Fig6a.eps', 'epsc')
close all;

%% Plot Fig 6b

N = 506;
Nalive = (N-1:-1:floor(0.02*N)-1)';
figure(2)
box on; hold on;

for chi = [17,31]
    infile = sprintf('Np_1_chi_%d_killTime_all.txt', chi);
    data_all1 = load(infile, '-ascii');
    infile = sprintf('Np_2_chi_%d_killTime_all.txt', chi);
    data_all2 = load(infile, '-ascii');
    
    kappa_all = data_all1./data_all2;
    
    for n_prey = 1:numel(Nalive)
        dat = kappa_all(n_prey, kappa_all(n_prey,:)>0);
        kappa_mean(n_prey,1) = mean(dat);
    end
    figure(2)
    plot(1-Nalive./N, kappa_mean)
end

figure(2)
plot([0 1],[2 2], '-.k')
xlabel('N_d'); ylabel('\kappa');
legend('\chi=17', '\chi=31', '\kappa=2');
title('Fig6b');
saveas(gcf, 'Fig6b.eps', 'epsc')