% CBCL

clear all; close all; clc;
cd('D:\MAVNMF\GitHub');
restoredefaultpath
addpath(genpath('D:\MAVNMF\GitHub'))
rng(1);
load('CBCL.mat');
[m,n] = size(X);
r = 49;
options.timemax = Inf;
options.maxiter = 1000;
options.model = 3;
% SPA init
options.display = 0; 

options.lambda = 0.89;
disp('Running min-vol NMF:'); 
[W_min,H_min,e_min,er1_min,er2_min] = minvolNMF(X,r,options);

% plot
affichage(W_min,7,19,19); 
exportgraphics(gcf,'plots/cbclmin.pdf','ContentType','vector');

options.model = 2;
options.lambda = 0.89;
disp('Running max-vol NMF:'); 
[W_max,H_max,e_max,er1_max,er2_max] = minvolNMF(X',r,options);
temp = H_max;
H_max = W_max';
W_max = temp';
% plot
affichage(W_max,7,19,19);
exportgraphics(gcf,'plots/cbclmax.pdf','ContentType','vector');

numZeros = nnz(W_min == 0)
numZeros = nnz(W_max == 0)

numZeros = nnz(H_min == 0)
numZeros = nnz(H_max == 0)

tol = 1e-10;
numZeros = nnz(abs(W_min) < tol)
numZeros = nnz(abs(W_max) < tol)

numZeros = nnz(abs(H_min) < tol)
numZeros = nnz(abs(H_max) < tol)

delta_default = 0.1;
vol_min_m = log(det ( W_min'*W_min + delta_default*eye(r) ));
vol_max_m = log(det ( W_max'*W_max + delta_default*eye(r) ));
vol_m = [vol_min_m, vol_max_m]';
round(vol_m, 3)

% min(min(H_min))
% min(min(H_max))
% 
% min(min(W_min))
% min(min(W_max))