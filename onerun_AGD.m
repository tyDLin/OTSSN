%%%

clear; 
clc; 

d = 5; 
sigma = 0.005; 
nsamples = [500]; 

for i = 1:1:length(nsamples)
    
    [X, Y, X_f, Y_f] = generate_data(nsamples(i), d); 
    data = kernel_make(X, Y, X_f, Y_f, sigma); 

    reg1 = 1/nsamples(i); 
    reg2 = 1/sqrt(nsamples(i)); 

    [~, ~, ~, res_time_AGD, res_norm_AGD] = AGD(data, reg1, reg2, 0); 
    [~, ~, ~, res_time_SSN, res_norm_SSN] = SSN(data, reg1, reg2, 0);

    %% plot the figures
    figure; 
    plot(res_time_AGD, res_norm_AGD, '-*', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(res_time_SSN, res_norm_SSN, '-d', 'LineWidth', 3, 'MarkerSize', 15);
    hold off
    legend('AGD (Gradient norm)', 'Specialized SSN (Residue norm)', 'Location', 'northeast', 'Orientation', 'vertical');

    grid on
    set(gca, 'YScale','log');
    set(gca, 'FontSize', 20);
    xlabel('Time (seconds)');
%    ylabel('Residue Norm');
    title(['(sample size, dimension) = (', num2str(nsamples(i)), ', ', num2str(d), ')']);

    path = sprintf('./plot_AGD_sample_%i', nsamples(i)); 
    saveas(gcf, path, 'epsc');
end