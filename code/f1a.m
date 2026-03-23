%% illustrate non-uniqueness
clear; close all; clc;
cd('D:\MAVNMF\GitHub');
restoredefaultpath
addpath(genpath('D:\MAVNMF\GitHub'))
%% synthetic data generation 
rng(1); % guarantee reproductibility

% Generate the matrix Wt of true basis vectors 
basis = [1, 0, 0;
      0, 1, 0;
      0, 0, 1];

Wt = rand(3, 3);

edge12 = points_on_line(basis(1,:), basis(2,:), 6/8);
edge21 = points_on_line(basis(1,:), basis(2,:), 1 - 7/8);

edge23 = points_on_line(basis(2,:), basis(3,:), 7/10);
edge32 = points_on_line(basis(2,:), basis(3,:), 1 - 8/10);

edge31 = points_on_line(basis(3,:), basis(1,:), 9/12);
edge13 = points_on_line(basis(3,:), basis(1,:), 1 - 10/12);

extra_points = [
    edge12;
    edge21;
    edge23;
    edge32;
    edge31;
    edge13
];

Wt = [
    Wt;
    extra_points
    ]; 

r=3;

% fix the number of data points
n = 50;

min_abundance = 0.2;
max_abundance = 0.8;
true_H = rand(r, n);
true_H = true_H ./ sum(true_H, 1);

currentMin = min(true_H, [], 'all');
currentMax = max(true_H, [], 'all');

true_H = min_abundance + (max_abundance - min_abundance) * (true_H - currentMin) / (currentMax - currentMin);

true_H = true_H ./ sum(true_H, 1);
Ht = true_H;
min(true_H, [], 2)
max(max(true_H))

X = Wt*Ht;

x1 = 0.5 * (2 * basis(2,:) + basis(3,:));
x2 = (sqrt(3)/2) * basis(3,:);
ternary_to_xy_Wt = [x1; x2]';
Wt_triangle = [ternary_to_xy_Wt; ternary_to_xy_Wt(1,:)];

x1 = 0.5 * (2 * Ht(2,:) + Ht(3,:));
x2 = (sqrt(3)/2) * Ht(3,:);
ternary_to_xy_X = [x1; x2];

%% min NMF
options=[];
options.display = 0;
options.maxiter = 1000;
options.model = 4;

options.lambda = 0.001;%0.001;S
disp('Running MVC-NMF:'); 
[W_min,H_min,e_min,er1_min,er2_min] = minvolNMF(X,r,options); 
delta_default = 0.1;
vol_min = log( det ( W_min'*W_min + delta_default*eye(r) ) );

fprintf('Min Vol Error ||Wmin-Wt||/||Wt|| = %2.2f%%.\n', 100*compareWs( Wt, W_min ) );


utransform = (W_min' * pinv(Wt'))';
x1 = 0.5 * (2 * utransform(2,:) + utransform(3,:));
x2 = (sqrt(3)/2) * utransform(3,:);
ternary_to_xy_W_min = [x1; x2]';

%% max NMF
options=[];
options.display = 0;
options.maxiter = 1000;
options.model = 5;

options.lambda = 0.02;%0.001;S
disp('Running MAV-NMF:'); 
[W_max,H_max,e_max,er1_max,er2_max] = minvolNMF(X',r,options); 
temp = H_max;
H_max = W_max';
W_max = temp';
delta_default = 0.1;
vol_max = log( det ( W_max'*W_max + delta_default*eye(r) ) );

fprintf('Max Vol Error ||Wmax-Wt||/||Wt|| = %2.2f%%.\n', 100*compareWs( Wt, W_max ) );


utransform = (W_max' * pinv(Wt'))';
x1 = 0.5 * (2 * utransform(2,:) + utransform(3,:));
x2 = (sqrt(3)/2) * utransform(3,:);
ternary_to_xy_W_max = [x1; x2]';

%% lambda min
W_min_triangle = [ternary_to_xy_W_min; ternary_to_xy_W_min(1,:)];

% lambda max
W_max_triangle = [ternary_to_xy_W_max; ternary_to_xy_W_max(1,:)];


Wt_triangle_2 = [    
    0.1    0.1
    0.6   0.8
    0.7931    0.0790
    0.1    0.1];

R = @(t)[cos(t) -sin(t); sin(t) cos(t)];

angles = 0.3;
scales = 0.9;
center_true = mean(Wt_triangle, 2);
Wt_triangle_3 = ( ...
    scales * R(angles) * (W_max_triangle - center_true)' ...
    + center_true' ...
)';
Wt_triangle_3(1, 1) = 0.08;
Wt_triangle_3(4, 1) = 0.08;

%% draw triangle
figure('Position', [100, 100, 800, 800]);
plot(Wt_triangle(:,1), Wt_triangle(:,2), 'b--', 'LineWidth', 5);
axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on
plot(Wt_triangle_2(:,1), Wt_triangle_2(:,2), 'b--', 'LineWidth', 5);
plot(Wt_triangle_3(:,1), Wt_triangle_3(:,2), 'b--', 'LineWidth', 5);

% data points
data_points_x = plot(ternary_to_xy_X(1,:), ternary_to_xy_X(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 5);

text(Wt_triangle(1,1)-0.24, Wt_triangle(1,2)-0.08, '$\textbf{M}_\textbf{1}\textbf{(:, 1)}$', 'Interpreter','latex','FontSize',35);
text(Wt_triangle(2,1)-0.1, Wt_triangle(2,2)-0.03, '$\textbf{M}_\textbf{1}\textbf{(:, 2)}$', 'Interpreter','latex','FontSize',35);
text(Wt_triangle(3,1),  Wt_triangle(3,2)+0.015, '$\textbf{M}_\textbf{1}\textbf{(:, 3)}$', 'Interpreter','latex','FontSize',35);

text(Wt_triangle_2(1,1)-0.2, Wt_triangle_2(1,2)+0.05, '$\textbf{M}_\textbf{2}\textbf{(:, 1)}$', 'Interpreter','latex','FontSize',35);
text(Wt_triangle_2(2,1),      Wt_triangle_2(2,2)+0.01, '$\textbf{M}_\textbf{2}\textbf{(:, 3)}$', 'Interpreter','latex','FontSize',35);
text(Wt_triangle_2(3,1)-0.02, Wt_triangle_2(3,2)+0.01, '$\textbf{M}_\textbf{2}\textbf{(:, 2)}$', 'Interpreter','latex','FontSize',35);


text(Wt_triangle_3(1,1)+0.05, Wt_triangle_3(1,2)-0.08, '$\textbf{M}_\textbf{3}\textbf{(:, 1)}$', 'Interpreter','latex','FontSize',35);
text(Wt_triangle_3(2,1)-0.2,      Wt_triangle_3(2,2)+0.03, '$\textbf{M}_\textbf{3}\textbf{(:, 2)}$', 'Interpreter','latex','FontSize',35);
text(Wt_triangle_3(3,1)-0.3, Wt_triangle_3(3,2)+0.05, '$\textbf{M}_\textbf{3}\textbf{(:, 3)}$', 'Interpreter','latex','FontSize',35);

exportgraphics(gcf,'plots/f1a.pdf','ContentType','vector');

vol = [vol_min, vol_max];
round(vol, 3)