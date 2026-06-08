%% illustrate min, max volume when data is highly mixed three dimensional data
clear; close all; clc;
cd('D:\MAVNMF\Revision 1');
restoredefaultpath
addpath(genpath('D:\MAVNMF\Revision 1'))
rng(1);

%%% load files
datasetFolder = 'bbcnewsdata\data\bbcnews';

files = dir(fullfile(datasetFolder,'**','*.txt'));

n = length(files);

documents = strings(n,1);
labels = strings(n,1);

for i = 1:n

    filename = fullfile(files(i).folder,files(i).name);

    documents(i) = string(fileread(filename));

    [~,label] = fileparts(files(i).folder);

    labels(i) = label;
end

labels = categorical(labels);

delta_default = 0.1;

cv = cvpartition(labels,'KFold',5);
acc_raw = NaN(cv.NumTestSets,1);
macroPrecision_raw = NaN(cv.NumTestSets,1);
macroRecall_raw    = NaN(cv.NumTestSets,1);
macroF1_raw        = NaN(cv.NumTestSets,1);


acc_sep = NaN(cv.NumTestSets,1);
macroPrecision_sep = NaN(cv.NumTestSets,1);
macroRecall_sep    = NaN(cv.NumTestSets,1);
macroF1_sep        = NaN(cv.NumTestSets,1);

acc_spa = NaN(cv.NumTestSets,1);
macroPrecision_spa = NaN(cv.NumTestSets,1);
macroRecall_spa    = NaN(cv.NumTestSets,1);
macroF1_spa        = NaN(cv.NumTestSets,1);

acc_min = NaN(cv.NumTestSets,1);
macroPrecision_min = NaN(cv.NumTestSets,1);
macroRecall_min    = NaN(cv.NumTestSets,1);
macroF1_min        = NaN(cv.NumTestSets,1);

acc_max = NaN(cv.NumTestSets,1);
macroPrecision_max = NaN(cv.NumTestSets,1);
macroRecall_max    = NaN(cv.NumTestSets,1);
macroF1_max        = NaN(cv.NumTestSets,1);

runtime_sep   = NaN(cv.NumTestSets,1);
runtime_spa  = NaN(cv.NumTestSets,1);
runtime_min  = NaN(cv.NumTestSets,1);
runtime_max   = NaN(cv.NumTestSets,1);

vol_sep   = NaN(cv.NumTestSets,1);
vol_spa  = NaN(cv.NumTestSets,1);
vol_min  = NaN(cv.NumTestSets,1);
vol_max   = NaN(cv.NumTestSets,1);
r = 100;

for fold = 1:cv.NumTestSets

    % Split data
    trainIdx = training(cv,fold);
    testIdx  = test(cv,fold);
   
    %%% Term documnet matrix
    documentsTrain = documents(trainIdx);
    documentsTest  = documents(testIdx);
    
    %% Preprocess training documents
    % documentsTrain = lower(string(documentsTrain));
    docsTrain = tokenizedDocument(documentsTrain);
    docsTrain = erasePunctuation(docsTrain);
    docsTrain = removeStopWords(docsTrain);
    docsTrain = normalizeWords(docsTrain,'Style','stem');
    
    %% Build training vocabulary
    bagTrain = bagOfWords(docsTrain);
    
    % Remove words appearing in fewer than 10 training documents
    bagTrain = removeInfrequentWords(bagTrain,10);
    
    %% Training document-term matrix
    Xtrain = full(bagTrain.Counts');
    
    %% Preprocess test documents
    % documentsTest  = lower(string(documentsTest));
    docsTest = tokenizedDocument(documentsTest);
    docsTest = erasePunctuation(docsTest);
    docsTest = removeStopWords(docsTest);
    docsTest = normalizeWords(docsTest,'Style','stem');
    
    %% Create test document-term matrix using training vocabulary
    Xtest = full(encode(bagTrain,docsTest))';    
 
   
    ytrain = labels(trainIdx);
    ytest  = labels(testIdx);

    %% RAW
    Mdl_raw = fitcecoc(Xtrain',ytrain);
    pred_raw = predict(Mdl_raw,Xtest');
    pred_raw = categorical(pred_raw);
    ytest = categorical(ytest);


    % Accuracy
    acc_raw(fold) = mean(pred_raw == ytest);

    fprintf('Fold %2d Accuracy_raw = %.4f\n', ...
        fold, acc_raw(fold));

    C_raw = confusionmat(ytest, pred_raw);
    
    numClasses_raw = size(C_raw,1);
    
    precision = NaN(numClasses_raw,1);
    recall    = NaN(numClasses_raw,1);
    f1        = NaN(numClasses_raw,1);

    for i = 1:numClasses_raw
        TP = C_raw(i,i);
        FP = sum(C_raw(:,i)) - TP;
        FN = sum(C_raw(i,:)) - TP;
    
        precision(i) = TP / (TP + FP);
        recall(i)    = TP / (TP + FN);
        f1(i)        = 2 * precision(i) .* recall(i) ./ (precision(i) + recall(i));
    end

    macroPrecision_raw(fold) = mean(precision);
    macroRecall_raw(fold)    = mean(recall);
    macroF1_raw(fold)        = mean(f1);
    fprintf('Fold %2d macroPrecision_raw = %.4f\n', ...
        fold, macroPrecision_raw(fold));
    fprintf('Fold %2d macroRecall_raw = %.4f\n', ...
        fold, macroRecall_raw(fold));
    fprintf('Fold %2d macroF1_raw = %.4f\n', ...
        fold, macroF1_raw(fold));


    %% SNPA
    optionssnpa.display = 0;
    optionssnpa.proj = 1;
    disp('Running SNPA NMF:'); 
    tStart = tic;
    [range,~] = SNPA(Xtrain,r,optionssnpa);
    W_sep = Xtrain(:,range);
    W_sep = W_sep./sum(W_sep, 1);
    H_train_sep = NaN(r,size(Xtrain,2));
    for j = 1:size(Xtrain,2)
        H_train_sep(:,j) = lsqnonneg(W_sep,Xtrain(:,j));
    end
    runtime_sep(fold) = toc(tStart);

    H_test_sep = NaN(r,size(Xtest,2));
    for j = 1:size(Xtest,2)
        H_test_sep(:,j) = lsqnonneg(W_sep,Xtest(:,j));
    end

    % ==========================
    % Train classifier
    % ==========================
    Mdl_sep = fitcecoc(H_train_sep',ytrain);
    pred_sep = predict(Mdl_sep,H_test_sep');
    pred_sep = categorical(pred_sep);
    ytest = categorical(ytest);


    acc_sep(fold) = mean(pred_sep == ytest);

    fprintf('Fold %2d Accuracy_sep = %.4f\n', ...
        fold, acc_sep(fold));

    C_sep = confusionmat(ytest, pred_sep);
    
    numClasses_sep = size(C_sep,1);
    
    precision = NaN(numClasses_sep,1);
    recall    = NaN(numClasses_sep,1);
    f1        = NaN(numClasses_sep,1);

    for i = 1:numClasses_sep
        TP = C_sep(i,i);
        FP = sum(C_sep(:,i)) - TP;
        FN = sum(C_sep(i,:)) - TP;
    
        precision(i) = TP / (TP + FP);
        recall(i)    = TP / (TP + FN);
        f1(i)        = 2 * precision(i) .* recall(i) ./ (precision(i) + recall(i));
    end

    macroPrecision_sep(fold) = mean(precision);
    macroRecall_sep(fold)    = mean(recall);
    macroF1_sep(fold)        = mean(f1);
    fprintf('Fold %2d macroPrecision_sep = %.4f\n', ...
        fold, macroPrecision_sep(fold));
    fprintf('Fold %2d macroRecall_sep = %.4f\n', ...
        fold, macroRecall_sep(fold));
    fprintf('Fold %2d macroF1_sep = %.4f\n', ...
        fold, macroF1_sep(fold));
    vol_sep(fold) = log(det ( W_sep'*W_sep + delta_default*eye(r) ));
    %% Sparse NMF
    optionssparse.FPGM = 1;
    optionssparse.timemax = Inf;
    optionssparse.maxiter = 1000;
    optionssparse.display = 0; 
    optionssparse.sW = 0.001;
    disp('Running Sparse NMF')
    tStart = tic;
    [W_spa,H_train_spa,es,ts] = sparseNMF(Xtrain,r,optionssparse);
    W_spa = W_spa./sum(W_spa, 1);
    runtime_spa(fold) = toc(tStart);

    H_train_spa = NaN(r,size(Xtrain,2));
    for j = 1:size(Xtrain,2)
        H_train_spa(:,j) = lsqnonneg(W_spa,Xtrain(:,j));
    end
    H_test_spa = NaN(r,size(Xtest,2));
    for j = 1:size(Xtest,2)
        H_test_spa(:,j) = lsqnonneg(W_spa,Xtest(:,j));
    end

    % ==========================
    % Train classifier
    % ==========================
    Mdl_spa = fitcecoc(H_train_spa',ytrain);
    pred_spa = predict(Mdl_spa,H_test_spa');
    pred_spa = categorical(pred_spa);

    acc_spa(fold) = mean(pred_spa == ytest);

    fprintf('Fold %2d Accuracy_spa = %.4f\n', ...
        fold, acc_spa(fold));

    C_spa = confusionmat(ytest, pred_spa);
    
    numClasses_spa = size(C_spa,1);
    
    precision = NaN(numClasses_spa,1);
    recall    = NaN(numClasses_spa,1);
    f1        = NaN(numClasses_spa,1);

    for i = 1:numClasses_spa
        TP = C_spa(i,i);
        FP = sum(C_spa(:,i)) - TP;
        FN = sum(C_spa(i,:)) - TP;
    
        precision(i) = TP / (TP + FP);
        recall(i)    = TP / (TP + FN);
        f1(i)        = 2 * precision(i) .* recall(i) ./ (precision(i) + recall(i));
    end

    macroPrecision_spa(fold) = mean(precision);
    macroRecall_spa(fold)    = mean(recall);
    macroF1_spa(fold)        = mean(f1);
    fprintf('Fold %2d macroPrecision_spa = %.4f\n', ...
        fold, macroPrecision_spa(fold));
    fprintf('Fold %2d macroRecall_spa = %.4f\n', ...
        fold, macroRecall_spa(fold));
    fprintf('Fold %2d macroF1_spa = %.4f\n', ...
        fold, macroF1_spa(fold));
    vol_spa(fold) = log(det ( W_spa'*W_spa + delta_default*eye(r) ));
    %% Min NMF
    options.timemax = Inf;
    options.maxiter = 1000;
    options.model = 3;
    options.display = 0; 
    
    options.lambda = 0.001;%0.8;
    disp('Running min-vol NMF:'); 
    tStart = tic;
    [W_min,H_train_min,e_min,er1_min,er2_min] = minvolNMF(Xtrain,r,options);
    runtime_min(fold) = toc(tStart);

    H_test_min = NaN(r,size(Xtest,2));

    for j = 1:size(Xtest,2)
        H_test_min(:,j) = lsqnonneg(W_min,Xtest(:,j));
    end

    Mdl_min = fitcecoc(H_train_min',ytrain);
    pred_min = predict(Mdl_min,H_test_min');
    pred_min = categorical(pred_min);

    acc_min(fold) = mean(pred_min == ytest);

    fprintf('Fold %2d Accuracy_min = %.4f\n', ...
        fold, acc_min(fold));

    C_min = confusionmat(ytest, pred_min);
    
    numClasses_min = size(C_min,1);
    
    precision = NaN(numClasses_min,1);
    recall    = NaN(numClasses_min,1);
    f1        = NaN(numClasses_min,1);

    for i = 1:numClasses_min
        TP = C_min(i,i);
        FP = sum(C_min(:,i)) - TP;
        FN = sum(C_min(i,:)) - TP;
    
        precision(i) = TP / (TP + FP);
        recall(i)    = TP / (TP + FN);
        f1(i)        = 2 * precision(i) .* recall(i) ./ (precision(i) + recall(i));
    end

    macroPrecision_min(fold) = mean(precision);
    macroRecall_min(fold)    = mean(recall);
    macroF1_min(fold)        = mean(f1);
    fprintf('Fold %2d macroPrecision_min = %.4f\n', ...
        fold, macroPrecision_min(fold));
    fprintf('Fold %2d macroRecall_min = %.4f\n', ...
        fold, macroRecall_min(fold));
    fprintf('Fold %2d macroF1_min = %.4f\n', ...
        fold, macroF1_min(fold));
    vol_min(fold) = log(det ( W_min'*W_min + delta_default*eye(r) ));
    %% Max NMF
    options.model = 2;
    options.lambda = 0.001;
    disp('Running max-vol NMF:'); 
    tStart = tic;
    [W_max,H_max,e_max,er1_max,er2_max] = minvolNMF(Xtrain',r,options);
    temp = H_max;
    H_max = W_max';
    W_max = temp';
    H_train_max = H_max;
    runtime_max(fold) = toc(tStart);
    
    H_test_max = NaN(r,size(Xtest,2));
    for j = 1:size(Xtest,2)
        H_test_max(:,j) = lsqnonneg(W_max,Xtest(:,j));
    end

    Mdl_max = fitcecoc(H_train_max',ytrain);
    pred_max = predict(Mdl_max,H_test_max');
    pred_max = categorical(pred_max);

    acc_max(fold) = mean(pred_max == ytest);

    fprintf('Fold %2d Accuracy_max = %.4f\n', ...
        fold, acc_max(fold));

    C_max = confusionmat(ytest, pred_max);
    
    numClasses_max = size(C_max,1);
    
    precision = NaN(numClasses_max,1);
    recall    = NaN(numClasses_max,1);
    f1        = NaN(numClasses_max,1);

    for i = 1:numClasses_max
        TP = C_max(i,i);
        FP = sum(C_max(:,i)) - TP;
        FN = sum(C_max(i,:)) - TP;
    
        precision(i) = TP / (TP + FP);
        recall(i)    = TP / (TP + FN);
        f1(i)        = 2 * precision(i) .* recall(i) ./ (precision(i) + recall(i));
    end

    macroPrecision_max(fold) = mean(precision);
    macroRecall_max(fold)    = mean(recall);
    macroF1_max(fold)        = mean(f1);
    fprintf('Fold %2d macroPrecision_max = %.4f\n', ...
        fold, macroPrecision_max(fold));
    fprintf('Fold %2d macroRecall_max = %.4f\n', ...
        fold, macroRecall_max(fold));
    fprintf('Fold %2d macroF1_max = %.4f\n', ...
        fold, macroF1_max(fold));
    vol_max(fold) = log(det ( W_max'*W_max + delta_default*eye(r) ));
end

fprintf('Accuracy_raw: %.3f ± %.3f\n', mean(acc_raw), std(acc_raw));
fprintf('Precision_raw: %.3f ± %.3f\n', mean(macroPrecision_raw), std(macroPrecision_raw));
fprintf('Recall_raw: %.3f ± %.3f\n', mean(macroRecall_raw), std(macroRecall_raw));
fprintf('F1-score_raw: %.3f ± %.3f\n', mean(macroF1_raw), std(macroF1_raw));


fprintf('Accuracy_sep: %.3f ± %.3f\n', mean(acc_sep), std(acc_sep));
fprintf('Precision_sep: %.3f ± %.3f\n', mean(macroPrecision_sep), std(macroPrecision_sep));
fprintf('Recall_sep: %.3f ± %.3f\n', mean(macroRecall_sep), std(macroRecall_sep));
fprintf('F1-score_sep: %.3f ± %.3f\n', mean(macroF1_sep), std(macroF1_sep));


fprintf('Accuracy_min: %.3f ± %.3f\n', mean(acc_min), std(acc_min));
fprintf('Precision_min: %.3f ± %.3f\n', mean(macroPrecision_min), std(macroPrecision_min));
fprintf('Recall_min: %.3f ± %.3f\n', mean(macroRecall_min), std(macroRecall_min));
fprintf('F1-score_min: %.3f ± %.3f\n', mean(macroF1_min), std(macroF1_min));

fprintf('Accuracy_spa: %.3f ± %.3f\n', mean(acc_spa), std(acc_spa));
fprintf('Precision_spa: %.3f ± %.3f\n', mean(macroPrecision_spa), std(macroPrecision_spa));
fprintf('Recall_spa: %.3f ± %.3f\n', mean(macroRecall_spa), std(macroRecall_spa));
fprintf('F1-score_spa: %.3f ± %.3f\n', mean(macroF1_spa), std(macroF1_spa));


fprintf('Accuracy_max: %.3f ± %.3f\n', mean(acc_max), std(acc_max));
fprintf('Precision_max: %.3f ± %.3f\n', mean(macroPrecision_max), std(macroPrecision_max));
fprintf('Recall_max: %.3f ± %.3f\n', mean(macroRecall_max), std(macroRecall_max));
fprintf('F1-score_max: %.3f ± %.3f\n', mean(macroF1_max), std(macroF1_max));



fprintf('Runtime_sep: %.3f ± %.3f\n', mean(runtime_sep), std(runtime_sep));
fprintf('Runtime_min: %.3f ± %.3f\n', mean(runtime_min), std(runtime_min));
fprintf('Runtime_spa: %.3f ± %.3f\n', mean(runtime_spa), std(runtime_spa));
fprintf('Runtime_max: %.3f ± %.3f\n', mean(runtime_max), std(runtime_max));

fprintf('vol_sep: %.3f ± %.3f\n', mean(vol_sep), std(vol_sep));
fprintf('vol_min: %.3f ± %.3f\n', mean(vol_min), std(vol_min));
fprintf('vol_spa: %.3f ± %.3f\n', mean(vol_spa), std(vol_spa));
fprintf('vol_max: %.3f ± %.3f\n', mean(vol_max), std(vol_max));

fprintf('%.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f \n', mean(acc_sep), std(acc_sep), mean(acc_min), std(acc_min), mean(acc_spa), std(acc_spa), mean(acc_max), std(acc_max));
fprintf('%.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f \n', mean(macroPrecision_sep), std(macroPrecision_sep), mean(macroPrecision_min), std(macroPrecision_min), mean(macroPrecision_spa), std(macroPrecision_spa),  mean(macroPrecision_max), std(macroPrecision_max));
fprintf('%.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f \n', mean(macroRecall_sep), std(macroRecall_sep), mean(macroRecall_min), std(macroRecall_min), mean(macroRecall_spa), std(macroRecall_spa),  mean(macroRecall_max), std(macroRecall_max));
fprintf('%.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f \n', mean(macroF1_sep), std(macroF1_sep), mean(macroF1_min), std(macroF1_min), mean(macroF1_spa), std(macroF1_spa),  mean(macroF1_max), std(macroF1_max));

fprintf('%.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f & %.3f ± %.3f \n', mean(vol_sep), std(vol_sep), mean(vol_min), std(vol_min), mean(vol_spa), std(vol_spa),  mean(vol_max), std(vol_max));

