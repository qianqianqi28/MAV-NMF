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
n = 500;

purity=[0.75 0.75 0.75];

alpha=[2; 2; 2];  % Parameter of the Dirichlet distribution
% generate randomly the abundance matrix true_H with a Dirichlet distribution
true_H = [sample_dirichlet(alpha,n)']; 
for j = 1 : n
    while (true_H(1,j) > purity(1)) || (true_H(2,j) > purity(2)) || (true_H(3,j) > purity(3))
        true_H(:,j) = sample_dirichlet(alpha,1)';
    end
end
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

options.lambda = 0.1;%0.001;S
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

%% plot
figure('Position', [100, 100, 800, 800]);

% true
w_true = plot(Wt_triangle(1:3,1), Wt_triangle(1:3,2), 'b+', 'MarkerSize', 20, 'LineWidth', 3);

axis equal;
axis off;
set(gca, 'FontSize', 20);
hold on


% lambda min
W_min_triangle = [ternary_to_xy_W_min; ternary_to_xy_W_min(1,:)];
w_min = plot(W_min_triangle(1:3,1), W_min_triangle(1:3,2), 'cx', 'MarkerSize', 20, 'LineWidth', 3);

% lambda max
W_max_triangle = [ternary_to_xy_W_max; ternary_to_xy_W_max(1,:)];
w_max = plot(W_max_triangle(1:3,1), W_max_triangle(1:3,2), 'rd', 'MarkerSize', 20, 'LineWidth', 3);


% draw triangle
plot(Wt_triangle(:,1), Wt_triangle(:,2), 'b-', 'LineWidth', 3);
plot(W_min_triangle(:,1), W_min_triangle(:,2), 'c--', 'LineWidth', 3);
plot(W_max_triangle(:,1), W_max_triangle(:,2), 'r--', 'LineWidth', 3);

% data points
data_points_x = plot(ternary_to_xy_X(1,:), ternary_to_xy_X(2,:), 'bo', 'MarkerSize', 8, 'LineWidth', 3);

legend([w_true, w_min, w_max, data_points_x], {'M - True', 'M - MVC-NMF', 'M - MAV-NMF', 'X - data points'})
exportgraphics(gcf,'plots/minnomaxdirechletzero.pdf','ContentType','vector');

vol = [vol_min, vol_max];
round(vol, 3)