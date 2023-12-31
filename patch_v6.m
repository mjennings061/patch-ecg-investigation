%% Short spaced lead patch investigation into orthogonal leads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% patch_v6.m
% Author: Michael Jennings
% Updates:
% V1 - First commit using leadAug_V12.m and horacek_STelevation_v16.m
% V2 - Plot top 5 leads as median BSPM contour
% V3 - Added orthogonal lead + visualision
% V4 - Evaluate patch design
% V5 - Augmenting leads using Kornreich dataset
% V6 - Use Kornreich dataset to test for MI
% 10/12/19 - MI detection on STAFF dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Import data
fileNames = importdata('fileNames.txt');
BalloonBSPMdata = importdata('BalloonBSPMdata.mat');
points = importdata('daltorso.pts');    %3d torso points file
face = importdata('daltorso.fac');  %triangulation data

%% Definitions
fs = 500; %sampling freq
jDelay = 40e-3/(1/fs);  %a 40ms delay based after the annotated jpoint

%% PART 1 - SELECT SSL
%% Select only patients who are positive responders (+specific vessels)
[~, patientNo] = dataHorResponders(fileNames);

%% Calculate STE for all possible patients and all lead combinations
rankedData = steRank(BalloonBSPMdata, patientNo);

%% Filter ranked leads below the MAX_LEAD_DISTANCE
MAX_LEAD_DISTANCE = 100; %maximum distance between leads
rankedDataSVL = filterLeadLength(points, face, rankedData, MAX_LEAD_DISTANCE);
node_p = rankedDataSVL(1,1);  %node+ is column 2
node_n = rankedDataSVL(1,2); %node- is column 1
clear rankedData;

%% Calculate mean and median ST elevation across 12-lead, VSLs and chosen SSL
steElev = calculateSTE(BalloonBSPMdata, patientNo, jDelay, node_p, node_n);

%% Plot a BSPM of the top lead(s)
noLeads = 1; %Number of top-performing leads to display
plotBSPM(BalloonBSPMdata, patientNo, points, face, rankedDataSVL, jDelay, noLeads);
hold on;
clear noLeads;

% Create orthogonal lead of the top lead
node_n2 = 212; %new orthogonal lead
node_p2 = 234;
plot3(points(212,1),points(212,2),points(212,3),'xw','markersize',10);
plot3(points(234,1),points(234,2),points(234,3),'xw','markersize',10);
hold off;

%% Plot the best SSL for the most ST-elevated patient
idxBest = steMaxIdx(BalloonBSPMdata, patientNo, jDelay, node_p, node_n); %ranked index of patients
plotSSL(BalloonBSPMdata, idxBest(2), node_p, node_n, node_p2, node_n2, fs); %plot the orthogonal SSLs
clear MAX_LEAD_DISTANCE fileNames steElev rankedDataSVL idxBest face points;
% clear BalloonBSPMdata patientNo;

%% PART 2 - GENERATE SSL
%% Transformation coefficients from 12-lead to patch using Kornreich data
% dataKornreich = [importdata('data_Normal_352.mat'),  ...
%                 importdata('data_MI_352.mat'),  ...
%                 importdata('data_LVH_352.mat')];

% Annotate the first (1,1) sample. 4=normal, 5=MI, 6=LVH
data_normal = importdata('data_Normal_352.mat');
for i = 1:length(data_normal)
    data_normal{i}(1,1) = 4;
end
data_MI = importdata('data_MI_352.mat');
for i = 1:length(data_MI)
    data_MI{i}(1,1) = 5;
end
data_LVH = importdata('data_LVH_352.mat');
for i = 1:length(data_LVH)
    data_LVH{i}(1,1) = 6;
end
dataKornreich = [data_normal, data_MI, data_LVH];

[dataTrain dataTest] = splitData(dataKornreich, 0.8);   %partition data 80/20%
c_steLead = getCoeffs(dataTrain, node_p, node_n);  %get transform coefficients for STE lead
c_orthLead = getCoeffs(dataTrain, node_p2, node_n2);  %get transform coefficients for orthogonal lead
clear i %data_normal data_MI data_LVH;

%% Transform lead and measure performance using generated coefficients
N = length(dataTest);
ssl_ste = cell(1,N);
ssl_orth = cell(1,N);
for i = 1:N
    [ssl_ste{i}, rmse_ste(i), cc_ste(i)] = getBspmSSL(dataTest{i}, c_steLead, node_p, node_n);  %create the ST-specific SSL
    [ssl_orth{i}, rmse_orth(i), cc_orth(i)] = getBspmSSL(dataTest{i}, c_orthLead, node_p2, node_n2); %create the orthogonal SSL
end

rmse_ste_median = median(rmse_ste); %median RMSE across the test dataset
cc_ste_median = median(abs(cc_ste)); %median correlation coefficient across the dataset

%plot real lead vs calculated lead
[~, pNum] = max(cc_ste); %choose the patient with the highest CC
% pNum = 60; %patient number in the testing data (dataTest)
figure()
subplot(2,1,1)
compareLeads(dataTest{pNum}(node_p+3,:) - dataTest{pNum}(node_n+3,:), ssl_ste{pNum}, fs); %plot ste-specific lead
subplot(2,1,2)
compareLeads(dataTest{pNum}(node_p2+3,:) - dataTest{pNum}(node_n2+3,:), ssl_orth{pNum}, fs); %plot orthogonal lead
clear pNum dataTrain dataTest dataKornreich i N cc_orth cc_ste ...
    rmse_orth rmse_ste;

%% PART 3 - CLASSIFICATION
%% Annotate data for classification
% 0=normal; 1=MI; 2=LVH

data_Horacek = BalloonBSPMdata(patientNo);
for i = 1:length(data_Horacek)
    if(mod(i,2) == 1) %if odd index (baseline ECG)
        data_Horacek{i}(1,1) = 0;
    else data_Horacek{i}(1,1) = 1; %if even index (PBI ECG)   
    end
end
for i = 1:length(data_normal)
    data_normal{i}(1,1) = 0; %Normals annotated as 4
end
for i = 1:length(data_MI)
    data_MI{i}(1,1) = 1;    % MI is 1
end
for i = 1:length(data_LVH)
    data_LVH{i}(1,1) = 2;   % LVH is 2
end
dataAnn = [data_Horacek data_normal data_MI data_LVH];
clear BalloonBSPMdata data_Horacek data_normal data_MI data_LVH;

%% Split annotated data
[dataTrain, annTrain, dataTest, annTest] = splitDataAnn(dataAnn,0.8); %split data 80/20 train/test

%% Generate the SSL and extract features
ssl_STE_train = {};
ssl_orth_train = {};
ssl_1_train = {};
ssl_2_train = {};
ssl_3_train = {};
ssl_4_train = {};
for i = 1:length(dataTrain)
    % annotations for each lead
    ssl_STE_train{i}(1,:) = dataTrain{i}(1,:);
    ssl_orth_train{i}(1,:) = dataTrain{i}(1,:);
    ssl_1_train{i}(1,:) = dataTrain{i}(1,:);
    ssl_2_train{i}(1,:) = dataTrain{i}(1,:);
    ssl_3_train{i}(1,:) = dataTrain{i}(1,:);
    ssl_4_train{i}(1,:) = dataTrain{i}(1,:);
    
    % get the SSLs, including additional bipolar leads from the patch
    ssl_STE_train{i}(2,:) = dataTrain{i}(node_p+3,:) - dataTrain{i}(node_n+3,:);  %get ste lead
    ssl_orth_train{i}(2,:) = dataTrain{i}(node_p2+3,:) - dataTrain{i}(node_n2+3,:);  %get orth lead
    ssl_1_train{i}(2,:) = dataTrain{i}(node_p2+3,:) - dataTrain{i}(node_n+3,:); 
    ssl_2_train{i}(2,:) = dataTrain{i}(node_n2+3,:) - dataTrain{i}(node_n+3,:);
    ssl_3_train{i}(2,:) = dataTrain{i}(node_p+3,:) - dataTrain{i}(node_n2+3,:); 
    ssl_4_train{i}(2,:) = dataTrain{i}(node_p2+3,:) - dataTrain{i}(node_p+3,:); 
    
%     for j = 1:length(dataTrain{i}(1,:)) %for all generated SSLs
%         if(dataTrain{i}(1,j) == 3) %find the j-point
%             ssl_STE_train_j(i) = j+jDelay; %get j point sample no
%             ssl_STE_train_STval(i) = ssl_STE_train{i}(ssl_STE_train_j(i)); %get the st value
%             if(i>1) %prevent zero indexing
%                 ssl_STE_train_slope(i) = diff(ssl_STE_train{i}(ssl_STE_train_j(i-1:i))); %get the slope
%             else ssl_STE_train_slope(i) = 0; %first value will be zero
%             end
%             ssl_STE_train_STmean(i) = mean(ssl_STE_train{i}(ssl_STE_train_j(i)-jDelay/2:ssl_STE_train_j(i)+jDelay/2)); %mean ST-value (40ms spread)
%             break;
%         end
%     end
end

% Test set
ssl_STE_test = {};
ssl_orth_test = {};
ssl_1_test = {};
ssl_2_test = {};
ssl_3_test = {};
ssl_4_test = {};
for i = 1:length(dataTest)
    % annotations for each lead
    ssl_STE_test{i}(1,:) = dataTest{i}(1,:);    %first row are annotations
    ssl_orth_test{i}(1,:) = dataTest{i}(1,:);
    ssl_1_test{i}(1,:) = dataTest{i}(1,:);
    ssl_2_test{i}(1,:) = dataTest{i}(1,:);
    ssl_3_test{i}(1,:) = dataTest{i}(1,:);
    ssl_4_test{i}(1,:) = dataTest{i}(1,:);
    
    ssl_STE_test{i}(2,:) = dataTest{i}(node_p+3,:) - dataTest{i}(node_n+3,:);  %get ste lead
    ssl_orth_test{i}(2,:) = dataTest{i}(node_p2+3,:) - dataTest{i}(node_n2+3,:);  %get orth lead
    ssl_1_test{i}(2,:) = dataTest{i}(node_p2+3,:) - dataTest{i}(node_n+3,:); 
    ssl_2_test{i}(2,:) = dataTest{i}(node_n2+3,:) - dataTest{i}(node_n+3,:); 
    ssl_3_test{i}(2,:) = dataTest{i}(node_p+3,:) - dataTest{i}(node_n2+3,:); 
    ssl_4_test{i}(2,:) = dataTest{i}(node_p2+3,:) - dataTest{i}(node_p+3,:); 
end

%combine features and annotations together for classification
% ssl_STE_train_packaged = [ssl_STE_train_STval; ssl_STE_train_slope; ssl_STE_train_STmean; annTrain];

% normalise features for classification
% for i = 1:3
%     ssl_STE_train_packagedNorm(i,:) = normalize(ssl_STE_train_packaged(i,:),'range');
% end
% ssl_STE_train_packagedNorm(4,:) = annTrain;