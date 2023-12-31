%% Split data into training/test
%Inputs:
% - data - Raw data 1xN cell array
% - trainRatio - Normalised amount of training data (0->1) e.g. 0.8 = 80% training data
%Outputs:
% - dataTrain - Training data
% - dataTest - Test data
function [dataTrain, dataTest] = splitData(data, trainRatio)
tic;
    %split data into the object c
    c = cvpartition(length(data),'HoldOut',trainRatio);
    idx = c.test'; %boolean indexing
    dataTrain = data(idx); %1 in the index is training data
    dataTest = data(~idx); % (logical indexing)
t = toc;
disp(['splitData ', num2str(t), ' seconds']);
end