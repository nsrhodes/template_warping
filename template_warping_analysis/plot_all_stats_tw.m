% Plot all histograms of summary stats
%% Housekeeping
clear all
clc
close all



%% Gaussian over hist

% % Generate 20 random data points from a normal distribution
% data = tstat_cor;  % You can modify this line to generate different data
% 
% % Create a histogram
% figure;
% histogram(data,3, 'Normalization', 'pdf');  % Normalize to get probability density
% set(gca,'FontSize',12,'FontName','Arial')
% % Hold the current figure to overlay the Gaussian fit
% hold on;
% 
% % Fit a Gaussian to the data
% mu = mean(data);        % Mean of the data
% sigma = std(data);      % Standard deviation of the data
% 
% % Create a range of x values for the Gaussian curve
% x = linspace(min(data)-1, max(data)+1, 100);
% gaussianFit = (1/(sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma).^2);
% 
% % Plot the Gaussian fit
% plot(x, gaussianFit, 'r-', 'LineWidth', 2);
% 
% % Add titles and labels
% title('Histogram with Gaussian Fit');
% xlabel('Data Points');
% ylabel('Probability Density');
% legend('Histogram', 'Gaussian Fit');
% 
% % Release the hold
% hold off;
% Specify the Excel file name
filename = 'R:\DRS-KidsOPM\Temp_warp_paper\PAPER\temp_warp_results_data.xlsx'; % Change this to your Excel file name

% Initialize figure for plotting

% Loop through each of the 8 plots
for i = 1:8
    % Load data from the Excel sheets
    dataSheet1 = xlsread(filename, 'With_coreg_error', sprintf('B%d:I%d', 1, 21)); % Read all 20 rows from each column
    dataSheet2 = xlsread(filename, 'No_coreg_error', sprintf('B%d:I%d', 1, 21)); % Read all 20 rows from each column

    % Combine the data from the current column for each sheet
    data1 = dataSheet1(:, i); % Take column pairs for each plot
    data2 = dataSheet2(:, i); % Take column pairs for each plot

    % Create a subplot for the current plot
    h = figure(i);
    set(h,'Position',[200 200 400 200])
    if i == 1
        % Plot histogram for sheet 1 in pale blue
        histogram(data1, 10, 'Normalization', 'pdf', 'FaceColor', [0.6 0.8 1], 'EdgeColor', 'none'); % Pale blue
        hold on;

        % Fit Gaussian for sheet 1
        mu1 = mean(data1);
        sigma1 = std(data1);
        x1 = linspace(min(data1)-1, max(data1)+1, 100);
        gaussianFit1 = (1/(sigma1 * sqrt(2 * pi))) * exp(-0.5 * ((x1 - mu1) / sigma1).^2);
        plot(x1, gaussianFit1, 'Color', [0.2 0.4 0.8], 'LineWidth', 2); % Darker blue

        % Add titles and labels
        xlabel('Distance between AAL peaks (mm)');
        ylabel('Probability Density');
        legend('Data', 'Gaussian Fit');
        set(gca,'FontSize',12,'FontName','Arial')
        xlim([0 max(data1)])
        

    else
    % Plot histogram for sheet 1
    histogram(data1, 10, 'Normalization', 'pdf', 'FaceColor', [1 0.6 0.4], 'EdgeColor', 'none'); % Pale orange
    hold on;
    set(gca,'FontSize',12,'FontName','Arial')
    % Fit Gaussian for sheet 1
    mu1 = mean(data1);
    sigma1 = std(data1);
    x1 = linspace(min(data1)-1, max(data1)+1, 100);
    gaussianFit1 = (1/(sigma1 * sqrt(2 * pi))) * exp(-0.5 * ((x1 - mu1) / sigma1).^2);
    plot(x1, gaussianFit1, 'Color', [1 0.4 0.2], 'LineWidth', 2); % Darker orange

    % Plot histogram for sheet 2
    histogram(data2, 10, 'Normalization', 'pdf', 'FaceColor', [0.6 0.4 1], 'EdgeColor', 'none'); % Pale purple

    % Fit Gaussian for sheet 2
    mu2 = mean(data2);
    sigma2 = std(data2);
    x2 = linspace(min(data2)-1, max(data2)+1, 100);
    gaussianFit2 = (1/(sigma2 * sqrt(2 * pi))) * exp(-0.5 * ((x2 - mu2) / sigma2).^2);
    plot(x2, gaussianFit2, 'Color', [0.4 0.2 1], 'LineWidth', 2); % Darker purple
    legend('','With coregistration error','','Without coregistration error');
    ylabel('Probability Density');
    % Add titles and labels
    if i == 2
        xlabel('Peak beta modulation separation (mm)');
        xlim([0 max(data1)])
    elseif i == 3
         xlabel('Peak beta modulation separation (mm)');
        xlim([0 max(data1)])
    elseif i == 4
        xlabel('Pseudo-T correlation')
        xlim([0 1])
        legend('','With coregistration error','','Without coregistration error','Location','northwest');
    elseif i == 5
       xlabel('Pseudo-T correlation')
        xlim([0 1])
        legend('','With coregistration error','','Without coregistration error','Location','northwest');

    elseif i == 6
        xlabel('TFS correlation')
        xlim([0 1]);
        legend('','With coregistration error','','Without coregistration error','Location','northwest');

     elseif i == 7
        xlabel('TFS correlation')
        xlim([0 1]);
        legend('','With coregistration error','','Without coregistration error','Location','northwest');

     elseif i == 8
        xlabel('Connectivity matrices correlation')
        xlim([0 1]);
        legend('','With coregistration error','','Without coregistration error','Location','northwest');

    end


    % Release hold for next subplot
    hold off;
    end
end
