clear variables; close all; clc;

L = 2.5;
prey_dia = 2*0.0136;

skip_vals = 10;
ignore_initialData = 201;

%% Fig 4a insets
box on; hold on;
for chi = [17, 24, 31]
    infile = sprintf('Fig4a_chi_%d_v2.txt', chi);
    prey_data = load(infile, '-ascii');
    prey_data(:,1) = prey_data(:,1).*sqrt(L/prey_dia);

    slope = zeros(size(prey_data,1),1);
    
    for i = skip_vals+1:size(prey_data,1)
        slope(i,1) = -(prey_data(i,2)-prey_data(i-skip_vals,2))/(prey_data(i,1)-prey_data(i-skip_vals,1));
    end
    
    slope(1:ignore_initialData) = [];
    
    T = prey_data(ignore_initialData+1:end,1);
    
    xlabel('time'); ylabel('rate of killing');
    
    slope_smooth = smooth(slope, 0.05, 'rloess');
    semilogx(T, slope_smooth)   
end
legend('\chi=17','\chi=24','\chi=31')
saveas(gcf,'Fig4a_insets.eps', 'epsc');
close all;

%% Fig 4b insets
box on; hold on;
for chi = [17,24,31]
    infile = sprintf('Fig4b_chi_%d_v2.txt', chi);
    prey_data = load(infile, '-ascii');
    prey_data(:,1) = prey_data(:,1).*sqrt(L/prey_dia);

    slope = zeros(size(prey_data,1),1);
    
    for i = skip_vals+1:size(prey_data,1)
        slope(i,1) = -(prey_data(i,2)-prey_data(i-skip_vals,2))/(prey_data(i,1)-prey_data(i-skip_vals,1));
    end
    
    slope(1:ignore_initialData) = [];
    
    T = prey_data(ignore_initialData+1:end,1);
    
    xlabel('time'); ylabel('rate of killing');
    
    slope_smooth = smooth(slope, 0.05, 'rloess');
    semilogx(T, slope_smooth)
end
legend('\chi=17','\chi=24','\chi=31')
saveas(gcf, 'Fig4b_insets.eps', 'epsc');