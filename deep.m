clear all; close all; clc;


%% Format SSLs (training and test)
data = importdata('ssl_ste_data.mat');
len_max = 500;  % max number of samples in a recording for zero padding
N_train = length(data.annTrain); %number of recordings
N_test = length(data.annTest); %number of recordings

% zero pad each recording and store to sig train
sig_train = [];
for subject = 1:N_train
    sig_train = [sig_train; data.ssl_STE_train{subject}(2,:) ...
            zeros(1,len_max - length(data.ssl_STE_train{subject}(2,:)))];
    if(data.annTrain(subject) == 2)
        annTrain(subject) = 0;
    else
        annTrain(subject) = data.annTrain(subject);
    end
end
% repeat for test data
sig_test = [];
for subject = 1:N_test
    sig_test = [sig_test; data.ssl_STE_test{subject}(2,:) ...
            zeros(1,len_max - length(data.ssl_STE_test{subject}(2,:)))];
    % change annotations with a 2 (LVH) to 0 (not MI)
    if(data.annTest(subject) == 2)
        annTest(subject) = 0;
    else
        annTest(subject) = data.annTest(subject);
    end
end

% prep the data and annotations into a separate file
XTrain = sig_train';
YTrain = categorical(annTrain);
XTest = sig_test';
YTest = categorical(annTest);
rng(1);
% sig=cellfun(@(x,y)[x;y],num2cell(Normal_v2th_random_test(1:2:end,:),2),num2cell(Normal_v2th_random_test(2:2:end,:),2),'UniformOutput',false);

%% Learning
layers = [ ...
    sequenceInputLayer(500)
     bilstmLayer(50,'OutputMode','sequence')
%     bilstmLayer(50,'OutputMode','last')
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer
    ];
rng(1);
options = trainingOptions('adam', ...
    'MaxEpochs',10, ...
    'MiniBatchSize', 50, ...
    'InitialLearnRate', 0.01, ...
    'GradientThreshold', 1, ...
    'plots','training-progress', ...
    'Verbose',false);
rng(1);
net = trainNetwork(XTrain,YTrain,layers,options);
rng(1);
[trainPred,scores ]= classify(net,XTest);
plotconfusion(YTest',trainPred','Training Accuracy');

[X,Y,T,AUC] = perfcurve(YTest',scores(:,2),'1');
plot(X,Y)
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by BLSTM')
