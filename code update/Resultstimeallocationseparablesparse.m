clear all; close all; clc;
cd('D:\MAVNMF\Revision 1');
restoredefaultpath
addpath(genpath('D:\MAVNMF\Revision 1'))
rng(1);

%% ------------------------------------------------------------
%% Build full matrix from the two tables (Male + Female)
%% ------------------------------------------------------------

% Define column names (activities)
activity_cols = {"paidwork","dom.work","caring","shopping","per.need","eating","sleeping","educat.","particip","soc.cont","goingout","sports","gardening","outside","tv-radio","reading","relaxing","other"};
% Define row names
row_info = {
  "M1275","M1280","M1285",...
  "M2575","M2580","M2585",...
  "M3575","M3580","M3585",...
  "M5075","M5080","M5085",...
  "M6575","M6580","M6585",...
  "F1275","F1280","F1285",...
  "F2575","F2580","F2585",...
  "F3575","F3580","F3585",...
  "F5075","F5080","F5085",...
  "F6575","F6580","F6585"};

%% ---------------------- MALE ----------------------
M1275 = [901,87,33,120,289,508,3737,1447,128,515,490,419,111,48,752,272,78,146];
M1280 = [769,157,28,138,294,528,3765,1455,101,505,396,436,102,41,815,256,56,240];
M1285 = [707,155,15,127,316,527,3744,1537,92,449,441,485,100,64,860,188,73,200];

M2575 = [2180,250,194,152,293,623,3380,124,129,609,382,269,173,69,700,366,64,124];
M2580 = [1992,269,206,157,316,649,3403,245,126,649,321,279,213,35,671,318,58,172];
M2585 = [1899,341,184,183,302,605,3397,208,143,599,391,271,231,67,812,243,57,148];

M3575 = [1901,249,99,173,351,660,3463,56,195,671,360,206,259,88,785,316,59,188];
M3580 = [2008,289,128,157,339,709,3445,90,156,593,240,280,238,45,804,343,44,170];
M3585 = [2093,331,136,185,332,650,3347,85,148,479,336,291,268,64,812,319,58,146];

M5075 = [1708,244,51,227,350,709,3560,18,122,603,237,209,256,116,921,468,79,203];
M5080 = [1357,337,54,221,364,744,3569,58,207,704,279,299,288,76,862,413,57,190];
M5085 = [1206,450,25,230,352,686,3533,46,272,554,264,316,309,112,1012,467,68,174];

M6575 = [176,617,124,273,365,763,3801,10,159,811,213,297,366,86,1161,477,157,223];
M6580 = [71,563,27,251,392,767,3871,43,192,671,220,403,312,117,1198,660,92,230];
M6585 = [95,636,38,264,383,707,3694,54,214,619,274,476,308,178,1233,578,104,225];

%% ---------------------- FEMALE ----------------------

F1275 = [723,494,135,208,359,536,3744,1163,125,592,364,348,90,32,594,292,73,208];
F1280 = [665,460,99,200,377,513,3777,1321,88,557,400,370,76,32,581,257,74,234];
F1285 = [564,397,86,223,387,495,3821,1436,80,527,396,352,86,41,702,207,63,214];

F2575 = [439,1342,635,347,311,593,3526,77,85,780,316,306,149,41,547,300,88,199];
F2580 = [471,1338,673,336,339,607,3532,115,115,776,270,352,131,32,497,275,54,167];
F2585 = [704,1147,651,336,337,572,3447,120,115,736,303,368,145,44,565,265,63,164];

F3575 = [299,1567,296,372,325,664,3567,104,133,694,225,335,198,37,622,356,76,207];
F3580 = [375,1605,309,347,346,633,3554,98,143,689,229,440,154,34,576,311,63,174];
F3585 = [412,1529,308,373,351,656,3444,68,196,699,277,453,170,45,582,307,59,153];

F5075 = [151,1600,83,376,367,601,3673,27,195,758,255,323,197,53,710,478,78,154];
F5080 = [153,1558,84,335,368,613,3701,30,179,810,212,504,190,41,644,390,53,212];
F5085 = [233,1487,82,352,385,595,3566,40,195,721,268,545,217,76,708,377,64,170];

F6575 = [11,1319,78,384,372,635,3849,6,108,929,219,297,169,37,888,485,63,230];
F6580 = [6,1409,154,292,453,665,3713,21,124,796,187,482,191,39,860,404,67,216];
F6585 = [19,1318,44,320,366,615,3675,23,139,749,202,579,169,52,1076,460,69,204];

%% ---------------------- Assemble matrix ----------------------

% Create a cell array containing all row vectors

all_rows = {
  M1275,M1280,M1285,...
  M2575,M2580,M2585,...
  M3575,M3580,M3585,...
  M5075,M5080,M5085,...
  M6575,M6580,M6585,...
  F1275,F1280,F1285,...
  F2575,F2580,F2585,...
  F3575,F3580,F3585,...
  F5075,F5080,F5085,...
  F6575,F6580,F6585};

X = cell2mat(all_rows');
row_sums = sum(X, 2);
60*24*7
X = X';
for i = 1:18
    fprintf(['\\text{%s}' repmat(' & %d', 1, 30) '\\\\\n'], ...
        activity_cols{i}, X(i,1:30));
end
W_dim1 = sum(X, 2)/sum(sum(X));
H_dim1 = ones(1, 30);
X = X ./ sum(X, 1);
basis = [1, 0, 0;
      0, 1, 0;
      0, 0, 1];
x1 = 0.5 * (2 * basis(2,:) + basis(3,:));
x2 = (sqrt(3)/2) * basis(3,:);
ternary_to_xy_Wt = [x1; x2]';
Wt_triangle = [ternary_to_xy_Wt; ternary_to_xy_Wt(1,:)];

r = 3; 

%% Separable NMF
optionssnpa.proj = 1;
disp('Running SNPA NMF:'); 
[range,~] = SNPA(X,r,optionssnpa);
W_sep = X(:,range);
W_sep = W_sep ./ sum(W_sep, 1);
optionsnnlsfpgm.proj = 3;
H_sep = nnls_FPGM(X,W_sep,optionsnnlsfpgm); 

%% sparse NMF
optionssparse.FPGM = 1;
optionssparse.timemax = Inf;
optionssparse.maxiter = 1000;
optionssparse.display = 0; 
optionssparse.sW = 0.62;
disp('Running Sparse NMF')
[W_spa,H_spa,es,ts] = sparseNMF(X,r,optionssparse);

W_spa = W_spa ./ sum(W_spa, 1);
H_spa = H_spa ./ sum(H_spa, 1);

%%
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
options.lambda = 0.01;
disp('Running max-vol NMF:'); 
[W_max,H_max,e_max,er1_max,er2_max] = minvolNMF(X',r,options);
temp = H_max;
H_max = W_max';
W_max = temp';
W_max = W_max(:, [3, 2, 1]);
H_max = H_max([3, 2, 1], :);


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

W_com = [W_dim1, W_sep, W_min, W_spa, W_max];
H_com = [H_dim1', H_sep', H_min', H_spa', H_max'];
%fprintf('\\[\\bm{W} = \\begin{array}{l|c|ccc|ccc}\n');
%fprintf(' & 1 & 1 & 2 & 3 & 1 & 2 & 3  \\\\\\hline\n');

for i = 1:18
    fprintf('\\text{%s} & %.4f & %.4f & %.4f& %.4f & %.4f & %.4f & %.4f& %.4f & %.4f & %.4f& %.4f & %.4f & %.4f\\\\\n', ...
        activity_cols{i}, W_com(i,1), W_com(i,2), W_com(i,3), W_com(i,4), W_com(i,5), W_com(i,6), W_com(i,7), W_com(i,8), W_com(i,9), W_com(i,10), W_com(i,11), W_com(i,12), W_com(i,13));
end
for i = 1:30
    fprintf('\\text{%s} & %.4f & %.4f & %.4f& %.4f & %.4f & %.4f & %.4f& %.4f & %.4f & %.4f& %.4f & %.4f & %.4f\\\\\n', ...
        row_info{i}, H_com(i,1), H_com(i,2), H_com(i,3), H_com(i,4), H_com(i,5), H_com(i,6), H_com(i,7), H_com(i,8), H_com(i,9), H_com(i,10), H_com(i,11), H_com(i,12), H_com(i,13));
end
%fprintf('\\end{array}\\]\n');

delta_default = 0.1;
vol_sep_m = log(det ( W_sep'*W_sep + delta_default*eye(r) ));
vol_min_m = log(det ( W_min'*W_min + delta_default*eye(r) ));
vol_spa_m = log(det ( W_spa'*W_spa + delta_default*eye(r) ));
vol_max_m = log(det ( W_max'*W_max + delta_default*eye(r) ));
vol_m = [vol_sep_m, vol_min_m, vol_spa_m, vol_max_m]';
round(vol_m, 3)

min(min(H_sep))
min(min(H_min))
min(min(H_spa))
min(min(H_max))

min(min(W_sep))
min(min(W_min))
min(min(W_spa))
min(min(W_max))