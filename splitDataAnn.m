%% Split data into training/test with annotation
%Inputs:
% - data_Horacek - Raw data 1xN cell array for horacek patients
% - data_LVH - LVH patients
% - data_MI - MI patients
% - data_Normal - Normal patients
% - trainRatio - Normalised amount of training data (0->1) e.g. 0.8 = 80% training data
%Outputs:
% - dataTrain - Training data
% - annTrain - Training annotations
% - dataTest - Test data
% - annTest - Test annotations
function [dataTrain, annTrain, dataTest, annTest] = splitDataAnn(data,trainRatio)
tic;
    %split data into the object c
    c = cvpartition(length(data),'HoldOut',trainRatio);
    idx = c.test'; %boolean indexing
    dataTrain = data(idx); %1 in the index is training data
    dataTest = data(~idx); % (logical indexing)
    for i = 1:length(dataTrain)
        if(dataTrain{i}(1,1) == 0)
            annTrain(i) = 0;
        elseif(dataTrain{i}(1,1) == 1)
            annTrain(i) = 1;
        elseif(dataTrain{i}(1,1) == 2)
            annTrain(i) = 2;
        end
    end
    for i = 1:length(dataTest)
        if(dataTest{i}(1,1) == 0)
            annTest(i) = 0;
        elseif(dataTest{i}(1,1) == 1)
            annTest(i) = 1;
        elseif(dataTest{i}(1,1) == 2)
            annTest(i) = 2;
        end
    end
    
t = toc;
disp(['splitDataAnn ', num2str(t), ' seconds']);
end