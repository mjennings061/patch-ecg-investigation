clc;clear;close all;
DEBUG = 1;   %to plot graphs
DEBUG_i = 3; %index of patient to plot in training data
fs = 500;

%% Import data
addpath('SSL test data') 
load('SSLs.mat');

%% Extract J-point data
j_n_train = zeros(1,length(annTrain));
j_STE_train = zeros(1,length(annTrain));
j_orth_train = zeros(1,length(annTrain));
j_1_train = zeros(1,length(annTrain));
j_2_train = zeros(1,length(annTrain));
j_3_train = zeros(1,length(annTrain));
j_4_train = zeros(1,length(annTrain));
for i = 1:length(annTrain)  %for each patient
    for j = 1:length(ssl_STE_train{i}(1,:)) %and each sample
        % Extract the j point annotation and amplitude
        if(ssl_STE_train{i}(1,j) == 3)
            j_n_train(i) = j;   %store the j point sample number
            j_STE_train(i) = ssl_STE_train{i}(2,j); % store the j-point amplitudes
            j_orth_train(i) = ssl_orth_train{i}(2,j);
            j_1_train(i) = ssl_1_train{i}(2,j);
            j_2_train(i) = ssl_2_train{i}(2,j);
            j_3_train(i) = ssl_3_train{i}(2,j);
            j_4_train(i) = ssl_4_train{i}(2,j);
            break;
        end
    end
end
%Transpose
j_n_train = j_n_train';     %store the j point sample number
j_STE_train = j_STE_train'; % store the j-point amplitudes
j_orth_train = j_orth_train';
j_1_train = j_1_train';
j_2_train = j_2_train';
j_3_train = j_3_train';
j_4_train = j_4_train';

% % J-points in a matrix (row1=jPoint sample no., row2=STE J-point ... row7 = SSL4 j-point)
jPoints = [j_n_train j_STE_train j_orth_train j_1_train j_2_train j_3_train j_4_train];
jPointTrain = zeros(length(annTrain),7);    %preallocation
jPointTrain(:,1:6) = jPoints(:,1:6);    %jpoints for each patient
jPointTrain(:,7) = annTrain;    %annotations (0=normal, 1=MI, 2=LVH)
% jPointTrain = jPointTrain'; %transpose for readability
jPointTrain(jPointTrain(:,7)==2,7) = 0; % change all 2 annotations to 0
% jPoints = array2table(jPointTrain,...
%     'VariableNames',{'STE' 'Orth' 'SSL1' 'SSL2' 'SSL3' 'SSL4' 'Class'});
% jPoints.Class = logical(jPoints.Class);
% for i = 1:length(jPoints)
%     if(jPointTrain(i,7) == 0 || jPointTrain(i,7) == 2)
%         jPointTrain(i,7) = 'F';
%     elseif(jPointTrain(i,7) == 1)
%         jPointTrain(i,7) = 'T';
%     end   
% end
csvwrite('j.csv',jPointTrain);

%% Transpose SSLS
annTrain = annTrain';
ssl_STE_train = ssl_STE_train';
ssl_orth_train = ssl_orth_train';
ssl_1_train = ssl_1_train';
ssl_2_train = ssl_2_train';
ssl_3_train = ssl_3_train';
ssl_4_train = ssl_4_train';

annTest = annTest';
ssl_STE_test = ssl_STE_test';
ssl_orth_test = ssl_orth_test';
ssl_1_test = ssl_1_test';
ssl_2_test = ssl_2_test';
ssl_3_test = ssl_3_test';
ssl_4_test = ssl_4_test';

%% Extract mutual information

%% PCA analysis

%% Data formatting for deep learning
%LVH is not MI, so set all 2 annotations to 0
annTrain(annTrain(:,1)==2,1) = 0; % change all 2 annotations to 0
annTest(annTest(:,1)==2,1) = 0; 
annTrain = categorical(annTrain); %change annotations to categorical
annTest = categorical(annTest);

%Get rid of first row (beat fiducial marking)
for i = 1:length(ssl_STE_train)
    ssl_STE_train{i}(1,:) = [];
    ssl_orth_train{i}(1,:) = [];
    ssl_1_train{i}(1,:) = [];
    ssl_2_train{i}(1,:) = [];
    ssl_3_train{i}(1,:) = [];
    ssl_4_train{i}(1,:) = [];
end

% [ssl_STE_train,annTrain] = segmentSignals(ssl_STE_train,annTrain);
%split signals according to classes
miX = ssl_STE_train(annTrain=='1');
miY = annTrain(annTrain=='1');
normalX = ssl_STE_train(annTrain=='0');
normalY = annTrain(annTrain=='0');

%% Model definition


%% Training


%% Review performance


%% Plot graphs
if(DEBUG)
    % Plot all SSLs
    figure;
    subplot(3,2,1)
    plot(ssl_STE_train{DEBUG_i});
    title('STE Sensitive'); grid on;
    subplot(3,2,2)
    plot(ssl_orth_train{DEBUG_i});
    title('Orthognal'); grid on;
    subplot(3,2,3)
    plot(ssl_1_train{DEBUG_i});
    title('SSL lead 1'); grid on;
    subplot(3,2,4)
    plot(ssl_2_train{DEBUG_i});
    title('SSL lead 2'); grid on;
    subplot(3,2,5)
    plot(ssl_3_train{DEBUG_i});
    title('SSL lead 3'); grid on;
    subplot(3,2,6)
    plot(ssl_4_train{DEBUG_i});
    title('SSL lead 4'); grid on;
    
    %Plot signal lengths
    figure;
    L = cellfun(@length,ssl_STE_train);
    h = histogram(L);
    title('Signal Lengths')
    xlabel('Length')
    ylabel('Count')
end

%% clear old variables
clear i fs 