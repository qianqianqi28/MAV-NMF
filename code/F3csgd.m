%% illustrate min, max volume when data is highly mixed three dimensional data
clear; close all; clc;
cd('D:\MAVNMF\GitHub');
restoredefaultpath
addpath(genpath('D:\MAVNMF\GitHub'))
rng(1);
T = readtable('AppendixBtestdatasets.xlsx', 'Sheet', 'test3data');
X = table2array(T(1:end, 2:end));
colnames = 1:size(X, 2);
rowSums = sum(X, 2);
sum(sum(X, 2))
X = X';
size(X)
sum(X, 1)

basis_temp = readtable('AppendixBtestdatasets.xlsx', 'Sheet', 'truebasis');
true_basis = table2array(basis_temp(1:end, 2:end));
true_basis = true_basis(1:25, 1:3);
true_basis = true_basis ./ sum(true_basis, 1);

collabel = [ ...
    "<0.430"
    "0.43-0.56"
    "0.56-0.74"
    "0.74-0.98"
    "0.98-1.29"
    "1.29-1.7"
    "1.7-2.24"
    "2.24-2.96"
    "2.96-3.91"
    "3.91-5.15"
    "5.15-6.8"
    "6.8-8.97"
    "8.97-11.84"
    "11.84-15.6"
    "15.6-20.6"
    "20.6-27.2"
    "27.2-35.9"
    "35.9-47.4"
    "47.4-62.5"
    "62.5-82.5"
    "82.5-109"
    "109-144"
    "144-190"
    "190-250"
    "250-330"
];



r = 3; 

options.timemax = Inf;
options.maxiter = 1000;

% SPA init
options.display = 0; 
%% min vol
options.model = 4;
options.lambda = 0.01;
disp('Running min-vol NMF:'); 
[W_min,H_min,e_min,er1_min,er2_min] = minvolNMF(X,r,options);

%% max vol
options.model = 5;
options.lambda = 0.0008;
disp('Running max-vol NMF:'); 
[W_max,H_max,e_max,er1_max,er2_max] = minvolNMF(X',r,options);
temp = H_max;
H_max = W_max';
W_max = temp';

figure('Position', [100, 100, 800, 800]);
h1 = plot(colnames, true_basis, 'b-+', 'LineWidth', 3);%

set(gca, 'FontSize', 20);
hold on
h2 = plot(colnames, W_min, 'c--x', 'LineWidth', 3);%
h4 = plot(colnames, W_max, 'r-.diamond', 'LineWidth', 3);%
for k = 1:3
    [ymax, idx] = max(true_basis(:,k));
    xloc = colnames(idx);
    yloc = ymax;

    text(xloc-2, yloc, sprintf('Basis %d', k), ...
        'FontSize', 26, 'FontWeight', 'bold', ...
        'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom');
end
lgd = legend([h1(1), h2(1), h4(1)], {'M - True', 'M - MVC-NMF', 'M - MAV-NMF'}, 'Location','best')

lgd.FontSize = 26;
xticks(colnames);
xticklabels(collabel);
xtickangle(45);  % optional: rotate labels for readability
set(gca, 'FontSize', 12, 'FontWeight', 'bold');   % increase tick label size
xlabel('Grain size (μm)', 'FontSize', 26);
ylabel('Volume content', 'FontSize', 26);
ylim([0 0.3])
xlim([0 26])
exportgraphics(gcf,'plots/sgdminplainmaxmtest2.pdf','ContentType','vector');


delta_default = 0.1;
vol_min_m = log(det ( W_min'*W_min + delta_default*eye(r) ));
vol_max_m = log(det ( W_max'*W_max + delta_default*eye(r) ));
vol_m = [vol_min_m, vol_max_m]';
round(vol_m, 3)

min(min(H_min))
min(min(H_max))


Tsample = readtable('AppendixBtestdatasets.xlsx', 'Sheet', 'test3abundance');

Xsample = table2array(Tsample(1:end-1, 2:end-2));
min(min(Xsample))
max(max(Xsample))
