clear variables; close all; clc;

%% Plot Fig 7a

for chi = [163, 229, 294]
    infile = sprintf('Fig7a_chi_%d.txt', chi);
    prey_data = load(infile, '-ascii');
   
    % Plot Fig.7a main
    figure(1)
    semilogx(prey_data(:,1), prey_data(:,2))
    hold on;
end

figure(1)
xlabel('\tau'); ylabel('N_l');
legend('\chi=163', '\chi=229', '\chi=294');
title('Fig7a');
saveas(gcf, 'Fig7a_main.eps', 'epsc')

%% Plot Fig 7b

for chi = [163, 229, 294]
    infile = sprintf('Fig7b_chi_%d.txt', chi);
    prey_data = load(infile, '-ascii');
    
    % Plot Fig.7a main
    figure(2)
    semilogx(prey_data(:,1), prey_data(:,2))
    hold on;
    
end

figure(2)
xlabel('\tau'); ylabel('N_l');
legend('\chi=163', '\chi=229', '\chi=294');
title('Fig7b');
saveas(gcf, 'Fig7b_main.eps', 'epsc')