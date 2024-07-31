clear variables; close all; clc;

skip_vals = 10;
ignore_initialData = 201;

%% Fig 7a insets
for chi = [163, 229, 294]
    infile = sprintf('Fig7a_chi_%d.txt', chi);
    prey_data = load(infile, '-ascii');
    
    slope = zeros(size(prey_data,1),1);
    
    for i = skip_vals+1:size(prey_data,1)
        slope(i,1) = -(prey_data(i,2)-prey_data(i-skip_vals,2))/(prey_data(i,1)-prey_data(i-skip_vals,1));
    end
    
    slope(1:ignore_initialData) = [];
    
    T = prey_data(ignore_initialData+1:end,1);
    
    xlabel('time'); ylabel('rate of killing');
    
    slope_smooth = smooth(slope, 0.05, 'rloess');
    semilogx(T, slope_smooth)
    
    saveas(gcf,'Fig7a_insets.eps', 'epsc');
end

%% Fig 7b insets
for chi = [163, 229, 294]
    infile = sprintf('Fig7b_chi_%d.txt', chi);
    prey_data = load(infile, '-ascii');
    
    slope = zeros(size(prey_data,1),1);
    
    for i = skip_vals+1:size(prey_data,1)
        slope(i,1) = -(prey_data(i,2)-prey_data(i-skip_vals,2))/(prey_data(i,1)-prey_data(i-skip_vals,1));
    end
    
    slope(1:ignore_initialData) = [];
    
    T = prey_data(ignore_initialData+1:end,1);
    
    xlabel('time'); ylabel('rate of killing');
    
    slope_smooth = smooth(slope, 0.05, 'rloess');
    semilogx(T, slope_smooth)
    saveas(gcf, 'Fig7b_insets.eps', 'epsc');
end