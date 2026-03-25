%% illustrate nonuniqueness
clear; close all; clc;
cd('D:\MAVNMF\GitHub');
restoredefaultpath
addpath(genpath('D:\MAVNMF\GitHub'))
%% synthetic data generation 
rng(1); % guarantee reproductibility

r=3;
n = 50;

% H is of r rows and n columns
true_H = 0.3 + (0.7 - 0.3) * rand(r,n); 
true_H = true_H./sum(true_H, 1);

% Transform H into simplex coordinates
x1 = 0.5 * (2 * true_H(2,:) + true_H(3,:));
x2 = (sqrt(3)/2) * true_H(3,:);
ternary_to_xy_X = [x1; x2];

% plot
figure('Position', [100, 100, 800, 800]);

data_points_x = plot(ternary_to_xy_X(1,:), ternary_to_xy_X(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 5);
axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on

true_M = [0 1 0.5;
    0.1 -0.1 0.9];

R = @(t)[cos(t) -sin(t); sin(t) cos(t)];

angles = [0.4, -0.4, 0.8];
scales = [0.5, 0.5, 0.5];

k = 1;
M = scales(k) * R(angles(k)) * (true_M - mean(true_M,2)) + mean(true_M,2);
plot([M(1,:) M(1,1)], ...
     [M(2,:) M(2,1)], ...
     '--','LineWidth',5,'Color','b')

text(M(1,1)-0.1, M(2,1)-0.02, '$\textbf{M}_\textbf{1}\textbf{(:, 1)}$', 'Interpreter','latex','FontSize',30);
text(M(1,2)-0.08,      M(2,2)-0.03, '$\textbf{M}_\textbf{1}\textbf{(:, 2)}$', 'Interpreter','latex','FontSize',30);
text(M(1,3), M(2,3)+0.02, '$\textbf{M}_\textbf{1}\textbf{(:, 3)}$', 'Interpreter','latex','FontSize',30);

k = 2;
M = scales(k) * R(angles(k)) * (true_M - mean(true_M,2)) + mean(true_M,2);
plot([M(1,:) M(1,1)], ...
     [M(2,:) M(2,1)], ...
     '--','LineWidth',5,'Color','b')

text(M(1,1)-0.08, M(2,1)+0.03, '$\textbf{M}_\textbf{2}\textbf{(:, 1)}$', 'Interpreter','latex','FontSize',30);
text(M(1,2),      M(2,2)+0.01, '$\textbf{M}_\textbf{2}\textbf{(:, 2)}$', 'Interpreter','latex','FontSize',30);
text(M(1,3)-0.02, M(2,3)+0.01, '$\textbf{M}_\textbf{2}\textbf{(:, 3)}$', 'Interpreter','latex','FontSize',30);


k = 3;
M = scales(k) * R(angles(k)) * (true_M - mean(true_M,2)) + mean(true_M,2);
plot([M(1,:) M(1,1)], ...
     [M(2,:) M(2,1)], ...
     '--','LineWidth',5,'Color','b')

text(M(1,1)-0.08, M(2,1), '$\textbf{M}_\textbf{3}\textbf{(:, 1)}$', 'Interpreter','latex','FontSize',30);
text(M(1,2)-0.05,      M(2,2)+0.03, '$\textbf{M}_\textbf{3}\textbf{(:, 2)}$', 'Interpreter','latex','FontSize',30);
text(M(1,3)-0.02, M(2,3)+0.01, '$\textbf{M}_\textbf{3}\textbf{(:, 3)}$', 'Interpreter','latex','FontSize',30);

exportgraphics(gcf,'plots/nonuniqueness.pdf','ContentType','vector');