% CBCL

clear all; close all; clc;
cd('D:\MAVNMF\Revision 1');
restoredefaultpath
addpath(genpath('D:\MAVNMF\Revision 1'))
rng(1);
load('CBCL.mat');
[m,n] = size(X);
r = 49;

%% Separable NMF
optionssnpa.proj = 1;
disp('Running SNPA NMF:'); 
[range,~] = SNPA(X,r,optionssnpa);
W_sep = X(:,range);
W_sep = W_sep ./ sum(W_sep, 1);
optionsnnlsfpgm.proj = 0;
H_sep = nnls_FPGM(X,W_sep,optionsnnlsfpgm); 
affichage(W_sep,7,19,19); 
exportgraphics(gcf,'plots/cbclsep.pdf','ContentType','vector');
%% sparse NMF 
optionssparse.FPGM = 1;
optionssparse.timemax = Inf;
optionssparse.maxiter = 1000;
optionssparse.display = 0; 
optionssparse.sW = 0.82;
disp('Running Sparse NMF')
[W_spa,H_spa,es,ts] = sparseNMF(X,r,optionssparse);
W_spa = W_spa./sum(W_spa, 1);
% Display the basis images
affichage(W_spa,7,19,19); %title('sparse NMF (0.85)');
exportgraphics(gcf,'plots/cbclspa.pdf','ContentType','vector');

%%
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

numZeros = nnz(W_sep == 0)
numZeros = nnz(W_min == 0)
numZeros = nnz(W_spa == 0)
numZeros = nnz(W_max == 0)

numZeros = nnz(H_sep == 0)
numZeros = nnz(H_min == 0)
numZeros = nnz(H_spa == 0)
numZeros = nnz(H_max == 0)

tol = 1e-10;
numZeros = nnz(abs(W_sep) < tol)
numZeros = nnz(abs(W_min) < tol)
numZeros = nnz(abs(W_spa) < tol)
numZeros = nnz(abs(W_max) < tol)

numZeros = nnz(abs(H_sep) < tol)
numZeros = nnz(abs(H_min) < tol)
numZeros = nnz(abs(H_spa) < tol)
numZeros = nnz(abs(H_max) < tol)

delta_default = 0.1;
vol_sep_m = log(det ( W_sep'*W_sep + delta_default*eye(r) ));
vol_min_m = log(det ( W_min'*W_min + delta_default*eye(r) ));
vol_spa_m = log(det ( W_spa'*W_spa + delta_default*eye(r) ));
vol_max_m = log(det ( W_max'*W_max + delta_default*eye(r) ));
vol_m = [vol_sep_m, vol_min_m, vol_spa_m, vol_max_m]';
round(vol_m, 3)

% min(min(H_min))
% min(min(H_max))
% 
% min(min(W_min))
% min(min(W_max))