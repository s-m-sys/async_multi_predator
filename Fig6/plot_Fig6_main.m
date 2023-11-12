clear variables; close all; clc;

%% Plot Fig 6a

for chi = [163, 229, 294]
    infile = sprintf('Fig6a_chi_%d.txt', chi);
    prey_data = load(infile, '-ascii');
   
    % Plot Fig.6a main
    figure(1)
    semilogx(prey_data(:,1), prey_data(:,2))
    hold on;
end

figure(1)
xlabel('\tau'); ylabel('N_l');
legend('\chi=163', '\chi=229', '\chi=294');
title('Fig6a');
saveas(gcf, 'Fig6a_main.eps', 'epsc')

%% Plot Fig 6b

for chi = [163, 229, 294]
    infile = sprintf('Fig6b_chi_%d.txt', chi);
    prey_data = load(infile, '-ascii');
    
    % Plot Fig.6a main
    figure(2)
    semilogx(prey_data(:,1), prey_data(:,2))
    hold on;
    
end

figure(2)
xlabel('\tau'); ylabel('N_l');
legend('\chi=163', '\chi=229', '\chi=294');
title('Fig6b');
saveas(gcf, 'Fig6b_main.eps', 'epsc')