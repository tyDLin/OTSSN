%%

clear; 
clc; 

ndims = [2, 5, 10]; 
ntrials = 10; 
nsamples = [50, 100, 200, 400]; 

N = length(nsamples); 
sigma = 0.005; 

result_IPM = zeros(N, 1);
result_SSN = zeros(N, 1); 

for k = 1:1:length(ndims)
    
    d = ndims(k); 
    
    for i = 1:1:length(nsamples)
    
        tmp_IPM = zeros(ntrials, 1); 
        tmp_SSN = zeros(ntrials, 1);  

        for j = 1:1:ntrials
            [X, Y, X_f, Y_f] = generate_data(nsamples(i), d); 
            data = kernel_make(X, Y, X_f, Y_f, sigma); 

            reg1 = 1/nsamples(i); 
            reg2 = 1/sqrt(nsamples(i)); 

            [~, ~, t_IPM] = IPM(data, reg1, reg2, 0); 
            [~, ~, t_SSN, ~, ~] = SSN(data, reg1, reg2, 0);

            tmp_IPM(j) = t_IPM; 
            tmp_SSN(j) = t_SSN;
            fprintf('%d\t', j);
        end

        result_IPM(i) = mean(tmp_IPM); 
        result_SSN(i) = mean(tmp_SSN);  
        fprintf('\n');
    end

    %% plot the figures
    figure; 
    plot(nsamples, result_IPM, '-*', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(nsamples, result_SSN, '-d', 'LineWidth', 3, 'MarkerSize', 15);
    hold off
    legend('Short-step IPM', 'Specialized SSN', 'Location', 'northwest', 'Orientation', 'vertical');

    grid on
    set(gca, 'FontSize', 20);
    xlabel('Number of filling samples');
    ylabel('Time (seconds)');
    title(['The dimension is ', num2str(d)]);

    path = sprintf('./time_%i', d); 
    saveas(gcf, path, 'epsc');
end