%% illustrate min, max volume when data is highly mixed three dimensional data
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

Wt = rand(12, 3);

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

min_abundance = 0.1;
max_abundance = 0.9;
true_H = rand(r, n);
true_H = true_H ./ sum(true_H, 1);

currentMin = min(true_H, [], 'all');
currentMax = max(true_H, [], 'all');

true_H = min_abundance + (max_abundance - min_abundance) * (true_H - currentMin) / (currentMax - currentMin);

true_H = true_H ./ sum(true_H, 1);
Ht = true_H;

X = Wt*Ht; 

x1 = 0.5 * (2 * basis(2,:) + basis(3,:));
x2 = (sqrt(3)/2) * basis(3,:);
ternary_to_xy_Wt = [x1; x2]';
Wt_triangle = [ternary_to_xy_Wt; ternary_to_xy_Wt(1,:)];



%% min NMF
options=[];
options.display = 0;
options.maxiter = 1000;
options.model = 4;

options.lambda = 0.018;
disp('Running min-vol NMF:'); 
[W_min,H_min,e_min,er1_min,er2_min] = minvolNMF(X,r,options); 

fprintf('Min Vol Error ||Wmin-Wt||/||Wt|| = %2.2f%%.\n', 100*compareWs( Wt, W_min ) );

x1 = 0.5 * (2 * H_min(2,:) + H_min(3,:));
x2 = (sqrt(3)/2) * H_min(3,:);
ternary_to_xy_Hmin = [x1; x2];

W_min_transpo = W_min';
W_min_transpo = W_min_transpo ./ sum(W_min_transpo, 1);
x1 = 0.5 * (2 * W_min_transpo(2,:) + W_min_transpo(3,:));
x2 = (sqrt(3)/2) * W_min_transpo(3,:);
ternary_to_xy_Wmin_transpo = [x1; x2];



%% plot min H
K = convhull(ternary_to_xy_Hmin(1,:), ternary_to_xy_Hmin(2,:));
cone_H_min = ternary_to_xy_Hmin(:, K);

figure('Position', [100, 100, 800, 800]);
% true
plot(Wt_triangle(:,1), Wt_triangle(:,2), 'k-', 'LineWidth', 5);

axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on

text(Wt_triangle(1, 1)-0.05,Wt_triangle(1, 2)-0.05,'$\textbf{e}_\textbf{1}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(2, 1)+0.02,Wt_triangle(2, 2)-0.05,'$\textbf{e}_\textbf{2}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(3, 1),Wt_triangle(3, 2)+0.05,'$\textbf{e}_\textbf{3}$','Interpreter','latex','FontSize',50)

cone_h = plot(cone_H_min(1,:), cone_H_min(2,:), 'b--', 'LineWidth', 5);

plot(ternary_to_xy_Hmin(1,:), ternary_to_xy_Hmin(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 5);

exportgraphics(gcf,'plots/minvolh.pdf','ContentType','vector');

%% plot min W
K = convhull(ternary_to_xy_Wmin_transpo(1,:), ternary_to_xy_Wmin_transpo(2,:));
cone_W_min = ternary_to_xy_Wmin_transpo(:, K);

figure('Position', [100, 100, 800, 800]);

plot(Wt_triangle(:,1), Wt_triangle(:,2), 'k-', 'LineWidth', 5);

axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on
text(Wt_triangle(1, 1)-0.05,Wt_triangle(1, 2)-0.05,'$\textbf{e}_\textbf{1}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(2, 1)+0.02,Wt_triangle(2, 2)-0.05,'$\textbf{e}_\textbf{2}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(3, 1),Wt_triangle(3, 2)+0.05,'$\textbf{e}_\textbf{3}$','Interpreter','latex','FontSize',50)

cone_w = plot(cone_W_min(1,:), cone_W_min(2,:), 'b--', 'LineWidth', 5);

plot(ternary_to_xy_Wmin_transpo(1,:), ternary_to_xy_Wmin_transpo(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 5);

exportgraphics(gcf,'plots/minvolm.pdf','ContentType','vector');

%% max NMF
options=[];
options.display = 0;
options.maxiter = 1000;
options.model = 5;

options.lambda = 0.018;
disp('Running max-vol NMF:'); 
[W_max,H_max,e_max,er1_max,er2_max] = minvolNMF(X',r,options); 
temp = H_max;
H_max = W_max';
W_max = temp';

fprintf('Max Vol Error ||Wmax-Wt||/||Wt|| = %2.2f%%.\n', 100*compareWs( Wt, W_max ) );

x1 = 0.5 * (2 * H_max(2,:) + H_max(3,:));
x2 = (sqrt(3)/2) * H_max(3,:);
ternary_to_xy_Hmax = [x1; x2];

W_max_transpo = W_max';
W_max_transpo = W_max_transpo ./ sum(W_max_transpo, 1);
x1 = 0.5 * (2 * W_max_transpo(2,:) + W_max_transpo(3,:));
x2 = (sqrt(3)/2) * W_max_transpo(3,:);
ternary_to_xy_Wmax_transpo = [x1; x2];

%% plot max H
K = convhull(ternary_to_xy_Hmax(1,:), ternary_to_xy_Hmax(2,:));
cone_H_max = ternary_to_xy_Hmax(:, K);

figure('Position', [100, 100, 800, 800]);
% true
plot(Wt_triangle(:,1), Wt_triangle(:,2), 'k-', 'LineWidth', 5);

axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on

text(Wt_triangle(1, 1)-0.05,Wt_triangle(1, 2)-0.05,'$\textbf{e}_\textbf{1}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(2, 1)+0.02,Wt_triangle(2, 2)-0.05,'$\textbf{e}_\textbf{2}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(3, 1),Wt_triangle(3, 2)+0.05,'$\textbf{e}_\textbf{3}$','Interpreter','latex','FontSize',50)

cone_h = plot(cone_H_max(1,:), cone_H_max(2,:), 'b--', 'LineWidth', 5);

plot(ternary_to_xy_Hmax(1,:), ternary_to_xy_Hmax(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 5);

exportgraphics(gcf,'plots/maxvolh.pdf','ContentType','vector');

%% plot max W
K = convhull(ternary_to_xy_Wmax_transpo(1,:), ternary_to_xy_Wmax_transpo(2,:));
cone_W_max = ternary_to_xy_Wmax_transpo(:, K);

figure('Position', [100, 100, 800, 800]);

plot(Wt_triangle(:,1), Wt_triangle(:,2), 'k-', 'LineWidth', 5);

axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on
text(Wt_triangle(1, 1)-0.05,Wt_triangle(1, 2)-0.05,'$\textbf{e}_\textbf{1}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(2, 1)+0.02,Wt_triangle(2, 2)-0.05,'$\textbf{e}_\textbf{2}$','Interpreter','latex','FontSize',50)
text(Wt_triangle(3, 1),Wt_triangle(3, 2)+0.05,'$\textbf{e}_\textbf{3}$','Interpreter','latex','FontSize',50)

cone_w = plot(cone_W_max(1,:), cone_W_max(2,:), 'b--', 'LineWidth', 5);

plot(ternary_to_xy_Wmax_transpo(1,:), ternary_to_xy_Wmax_transpo(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 5);

exportgraphics(gcf,'plots/maxvolm.pdf','ContentType','vector');


vol_min_m = det ( W_min'*W_min );
vol_max_m = det ( W_max'*W_max );

vol_min_h = det ( H_min*H_min' );
vol_max_h = det ( H_max*H_max' );

vol_m = [vol_min_m, vol_max_m];
round(vol_m, 3)

vol_h = [vol_min_h, vol_max_h];
round(vol_h, 3)

vol_m.*vol_h

delta_default = 0.1;

vol_min_m = log(det ( W_min'*W_min  + delta_default*eye(r)));
vol_max_m = log(det ( W_max'*W_max  + delta_default*eye(r)));

vol_min_h = log(det ( H_min*H_min' + delta_default*eye(r)));
vol_max_h = log(det ( H_max*H_max' + delta_default*eye(r)));

vol_m = [vol_min_m, vol_max_m];
round(vol_m, 3)

vol_h = [vol_min_h, vol_max_h];
round(vol_h, 3)

vol_m.*vol_h