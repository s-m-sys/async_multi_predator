clear variables; close all; clc;

%% Plot Fig 4a

for chi = [17,24,31]
    infile = sprintf('Fig4a_chi_%d_v2.txt', chi);
    prey_data = load(infile, '-ascii');

    % Plot Fig.4a main
    figure(1)
    semilogx(prey_data(:,1), prey_data(:,2))
    hold on;
end

figure(1)
xlabel('\tau'); ylabel('N_l');
legend('\chi=17', '\chi=24', '\chi=31');
title('Fig4a');
xlim([1 2e6])
saveas(gcf, 'Fig4a_main.eps', 'epsc')

%% Plot Fig 4b

for chi = [17,24,31]
    infile = sprintf('Fig4b_chi_%d_v2.txt', chi);
    prey_data = load(infile, '-ascii');

    % Plot Fig.4b main
    figure(2)
    semilogx(prey_data(:,1), prey_data(:,2))
    hold on;
    
end

figure(2)
xlabel('\tau'); ylabel('N_l');
legend('\chi=17', '\chi=24', '\chi=31');
title('Fig4b');
xlim([1 2e6])
saveas(gcf, 'Fig4b_main.eps', 'epsc')