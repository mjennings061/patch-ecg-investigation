% Use horacek data extract 12-lead ECG info in beats_xxxx_Ho. Test the STEMI detector using this 
clc;clear;close all;
%% Select only patients who are positive responders (+specific vessels)
BalloonBSPMdata = importdata('BalloonBSPMdata.mat');
fileNames = importdata('fileNames.txt');
j=1;

% make an array of all patients and filter non-responders
for i=1:length(fileNames) 
    if (contains(fileNames(i),'_Y_')) %Only look at responders
        if(contains(fileNames(i), '_B_'))
            if(contains(fileNames(i+1),'_P_'))
                patientNo(j,1) = i; %patient number is stored for future reference
                patients(j,1) = fileNames(i);
                j=j+1;
                patientNo(j,1) = i+1;
                patients(j,1) = fileNames(i+1);
                if(contains(fileNames(i),'_F_'))
                    maleFemale(i:i+1) = 1;
                else
                    maleFemale(i:i+1) = 0;
                end
                j=j+1;
                i=i+2;
            end
        end
    end
end

% for i=1:length(fileNames) 
%     if (contains(fileNames(i),'_Y_') && contains(fileNames(i),'RCA')) %Only look at responders
%         if(contains(fileNames(i), '_B_'))
%             if(contains(fileNames(i+1),'_P_'))
%                 patientNo(j,1) = i; %patient number is stored for future reference
%                 patients(j,1) = fileNames(i);
%                 j=j+1;
%                 patientNo(j,1) = i+1;
%                 patients(j,1) = fileNames(i+1);
%                 j=j+1;
%                 i=i+2;
%             end
%         end
%     end
% end

%% Extract 9-lead ECG for each patient (V1-V6; I-III)
base_MF = ones(1,length(patientNo/2))*2;
for i = 1:2:length(patientNo)-1 %separate baseline and PBI patients
    PATIENT_NO = patientNo(i,1);
    base{i} = BalloonBSPMdata{1,PATIENT_NO}(:,:); %Baseline data
    infl{i} = BalloonBSPMdata{1,PATIENT_NO+1}(:,:);%Peak balloon data
    base_file(i) = patients(i);
    infl_file(i) = patients(i+1);
    base_MF(i) = maleFemale(PATIENT_NO); %male/female annotations
end
% infl_MF = base_MF; %sex is the same for both base and infl patient

%delete empty cells and filenames
base = base(~cellfun('isempty',base));
infl = infl(~cellfun('isempty',infl));
base_file = base_file(~cellfun('isempty',base_file));
infl_file = infl_file(~cellfun('isempty',infl_file));
base_MF(base_MF(:)==2) = [];
ssl_ste_base = cell(1,length(base));    % preallocation
ssl_ste_infl = cell(1,length(infl));    % preallocation
ssl_orth_base = cell(1,length(base));    % preallocation
ssl_orth_infl = cell(1,length(infl));    % preallocation

%calculate the 9 leads
for i = 1:length(base)
    beats_base_H{i}(:,1) = base{i}(169+3,:)'; %V1
    beats_base_H{i}(:,2) = base{i}(171+3,:)'; %V2
    beats_base_H{i}(:,3) = (base{i}(192+3,:)+base{i}(193+3,:))/2'; %V3
    beats_base_H{i}(:,4) = base{i}(216+3,:)'; %V4
    beats_base_H{i}(:,5) = (base{i}(217+3,:)+2*base{i}(218+3,:))/3'; %V5
    beats_base_H{i}(:,6) = base{i}(219+3,:)'; %V6
    beats_base_H{i}(:,7) = base{i}(2,:)'; %I
    beats_base_H{i}(:,8) = base{i}(3,:)'; %II
    beats_base_H{i}(:,9) = (base{i}(3,:)-base{i}(2,:)*2)'; %III
    
    beats_infl_H{i}(:,1) = infl{i}(169+3,:)'; %V1
    beats_infl_H{i}(:,2) = infl{i}(171+3,:)'; %V2
    beats_infl_H{i}(:,3) = (infl{i}(192+3,:)+infl{i}(193+3,:))/2'; %V3
    beats_infl_H{i}(:,4) = infl{i}(216+3,:)'; %V4
    beats_infl_H{i}(:,5) = (infl{i}(217+3,:)+2*infl{i}(218+3,:))/3'; %V5
    beats_infl_H{i}(:,6) = infl{i}(219+3,:)'; %V6
    beats_infl_H{i}(:,7) = infl{i}(2,:)'; %I
    beats_infl_H{i}(:,8) = infl{i}(3,:)'; %II
    beats_infl_H{i}(:,9) = (infl{i}(3,:)-infl{i}(2,:))*2'; %III
    
    %get the SSLs
    ssl_ste_base{i}(:,1) = base{i}(254+3,:) - base{i}(173+3,:);
    ssl_ste_infl{i}(:,1) = infl{i}(254+3,:) - infl{i}(173+3,:);
    ssl_orth_base{i}(:,1) = base{i}(234+3,:) - base{i}(212+3,:);
    ssl_orth_infl{i}(:,1) = infl{i}(234+3,:) - infl{i}(212+3,:);
    
    %baseline - get the jpoint and double it for oversampling
    for j = 1:length(base{i}(1,:))
        if(base{i}(1,j) == 3)
            jPoint_base(i) = j*2;
            break;
        end
    end
    %PBI - get the jpoint and double it for oversampling
    for j = 1:length(infl{i}(1,:))
        if(infl{i}(1,j) == 3)
            jPoint_infl(i) = j*2;
            break;
        end
    end
    
    %oversample via interpolation
    for j = 1:9
        beats_base_Ho{i}(:,j) = interp(beats_base_H{i}(:,j),2); %double sample rate to 1000Hz
        beats_infl_Ho{i}(:,j) = interp(beats_infl_H{i}(:,j),2);
    end
end

%% Get jpoints for SSLs
%col1 = st-sensitive lead; col2 = orthogonal lead; col3 = annotation (1 = STEMI)
ssl_jPoints = zeros(length(base)*2,3);
for i = 1:length(base)  % get j-points for baseline and inflation patients
    ssl_jPoints(i,1) = ssl_ste_base{i}(jPoint_base(i)/2);   % STE sensitive lead
    ssl_jPoints(i,2) = ssl_orth_base{i}(jPoint_base(i)/2);  % orthogonal lead
    ssl_jPoints(i,3) = 0;   % baseline (no STEMI)
    ssl_jPoints(length(base)+i,1) = ssl_ste_infl{i}(jPoint_infl(i)/2);   % STE sensitive lead
    ssl_jPoints(length(base)+i,2) = ssl_orth_infl{i}(jPoint_infl(i)/2);  % orthogonal lead
    ssl_jPoints(length(base)+i,3) = 1;   % STEMI
end

%% Classify stemi
age = 60;   %generic age
jPoint_base = jPoint_base+40;   %move jpoint by 40ms
jPoint_infl = jPoint_infl+40;
tp=0;fp=0;fn=0;tn=0;
for i = 1:length(base)  %for all patients
    beats_base_Ho{i} = beats_base_Ho{i}/1000;   %values are in uV, divide to mV
    beats_infl_Ho{i} = beats_infl_Ho{i}/1000;
    [result_base(i), location_base(i)] = detectSTEMI(beats_base_Ho{i}, jPoint_base(i), maleFemale(1), age); %classify STEMI
    [result_infl(i), location_infl(i)] = detectSTEMI(beats_infl_Ho{i}, jPoint_infl(i), maleFemale(1), age);
    
    %count tp,fp,tp,tn
    if(result_base(i) == 0)
        tn = tn+1;
    else
        fp = fp+1;
    end
    if(result_infl(i) == 0)
        fn = fn+1;
    else
        tp = tp+1;
    end
end

%% Calculate performance
sens = tp/(tp+fn)*100;  %sensitivity
spec = tn/(tn+fp)*100; %specificity
ppv = tp/(tp+fp)*100;%positive predictive value