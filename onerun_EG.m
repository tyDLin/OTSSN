% %%

clear; 
clc; 

d = 5; 
sigma = 0.005; 
% nsamples = [50, 100, 200]; 
nsamples = [500, 1000, 2000]; 

for i = 1:1:length(nsamples)
    
    [X, Y, X_f, Y_f] = generate_data(nsamples(i), d); 
    data = kernel_make(X, Y, X_f, Y_f, sigma); 

    reg1 = 1/nsamples(i); 
    reg2 = 1/sqrt(nsamples(i)); 

    [~, ~, ~, res_time_EG, res_norm_EG] = EG(data, reg1, reg2, 0); 
    [~, ~, ~, res_time_SSN, res_norm_SSN] = SSN(data, reg1, reg2, 0);

    %% plot the figures
    figure; 
    plot(res_time_EG, res_norm_EG, '-*', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(res_time_SSN, res_norm_SSN, '-d', 'LineWidth', 3, 'MarkerSize', 15);
    hold off
    legend('EG', 'Specialized SSN', 'Location', 'northwest', 'Orientation', 'horizontal');

    grid on
    set(gca, 'YScale','log');
    set(gca, 'FontSize', 20);
    xlabel('Time (seconds)');
    ylabel('Residue Norm');
    title(['The sample size is ', num2str(nsamples(i))]);

    path = sprintf('./plot_EG_samples_%i', nsamples(i)); 
    saveas(gcf, path, 'epsc');
end