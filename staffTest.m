warning('off','signal:filtfilt:ParseSOS');

% Split STAFF dataset into individual beats
clc;clear;close all;
DEBUG = 1; %plot debug patients
fs=1000;

%% Import STAFF annotations
opts = spreadsheetImportOptions("NumVariables", 9);
% Specify sheet and range
opts.Sheet = "Sheet";
opts.DataRange = "A2:I531";
% Specify column names and types
opts.VariableNames = ["filename", "patient", "age", "gender", "prior_mi", "artery", "inflation_start", "inflation_duration", "inflation_after"];
opts.VariableTypes = ["string", "double", "double", "categorical", "categorical", "categorical", "string", "string", "string"];
% Specify variable properties
opts = setvaropts(opts, ["filename", "inflation_start", "inflation_duration", "inflation_after"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["filename", "gender", "prior_mi", "artery", "inflation_start", "inflation_duration", "inflation_after"], "EmptyFieldRule", "auto");
% Import the data
ann = readtable("STAFF\ann.xlsx", opts, "UseExcel", false);
% Clear temporary variables
clear opts

%% Import files and separate annotations
path = 'STAFF/'; %folder with data
files = strcat(path,ann.filename,'.mat'); %extract from annotation file
j=1;k=1;
for i = 1:length(files) %for all files
    if(ann.prior_mi(i) == "no") %exclude prior MI
        if(ann.inflation_start(i) == "")    %exclude inflation ECGs
            filesBaseline(j) = ann.filename(i); %record baseline filename
            j = j+1;
        else
            ecg = importdata(files(i)); % import ECG file for that patient
            filesInflation(k) = ann.filename(i);    %record MI filename
            inflationStart_n(k) = str2num(ann.inflation_start(i))*fs;   %note inflation sample num
            if(ann.inflation_after(i) == "0")
                inflationEnd_n(k) = length(ecg(:,1)); %note the last sample no
            else
                %end of inflation 
                inflationEnd_n(k) = inflationStart_n(k) + str2num(ann.inflation_duration(i))*fs;
            end
            k = k+1;
        end
    end
end
inflationStartEnd_n = [inflationStart_n; inflationEnd_n]; %combine into start and end sample numbers
filesBaseline = strcat(path,filesBaseline,'.mat'); 
filesInflation = strcat(path,filesInflation,'.mat');
clear i j k path inflationEnd_n inflationStart_n ecg;

%% Baseline patients -> Extract beats
PATIENT_DEBUG = 2; %patient data to plot
disp('**** Processing baseline average complexes ****');
for i = 1:length(filesBaseline)
    %     i = 10; %override i, remove when testing on multiple patients
    %% Exclusion criteria
%     if(i==11 || i==17 || i==51 || i==53 || i==60 || i==70 || i==64 || ...
%             i==77 || i==134 || i==143 || i==144 || i==149 || i==150 || ...
%         i==151 || i==206 || i==241 || i==244)
%         beats_baseline{i} = [];
%         disp(['WARNING: Patient ', num2str(i), ' ignored due to short average beats']);
%         continue;
%     elseif(i==67 || i==74)
%         beats_baseline{i} = [];
%         disp(['WARNING: Patient ', num2str(i), ' ignored due to low SNR']);
%         continue;
%     end
    
    %% Import signal
    ecgsig = importdata(filesBaseline(i));  %import all leads for annotation (V1-6; I-III; 9 in total)]
    if(find(~isfinite(ecgsig))) %check if the data is finite
        beats_baseline{i} = [];
%         winStart_base(i) = 0;
%         winEnd_base(i) = 0;
        disp(['WARNING: Patient ', num2str(i), ' ignored due to non-finite elements']);
        continue;
    end
    
    %% Average complexes (plot graphs for one patient)
    if((i == PATIENT_DEBUG) && DEBUG)
%         [beats_baseline{i},winStart_base(i),winEnd_base(i)] = averageBeat(ecgsig, MIN_GAP, MIN_BEAT_LENGTH,1);
        beats_baseline{i} = averageBeat(ecgsig,1);
    else
%         [beats_baseline{i},winStart_base(i),winEnd_base(i)] = averageBeat(ecgsig, MIN_GAP, MIN_BEAT_LENGTH,0);
        beats_baseline{i} = averageBeat(ecgsig,0);
    end
    disp(['Patient ', num2str(i), ' completed'])
end
clear ecgsigAll ecgsig i j h 
disp('**** Baseline complexes complete ****');
disp('');

%% Inflation patients -> Extract beats
PATIENT_DEBUG = 3; %patient data to plot
MIN_GAP = 500; %minimum number of samples between beats
MIN_BEAT_LENGTH = 500; %minimum length of beat without error

disp('**** Processing inflation average complexes ****');
for i = 1:length(filesInflation)
    %% Exclusion criteria
%     i = 10; %override i, remove when testing on multiple patients
%     if(i==11 || i==55 || i==56 || i==60 || i==88)
%         beats_inflation{i} = [];
%         disp(['WARNING: Patient ', num2str(i), ' ignored due to short average beats']);
%         continue;
%     elseif(i==25 || i==53 || i==83)
%         beats_inflation{i} = [];
%         disp(['WARNING: Patient ', num2str(i), ' ignored due to low SNR']);
%         continue;
%     end
    
    %% ECG signal formatting
    ecgsig = importdata(filesInflation(i));  %import all leads for annotation (V1-6; I-III; 9 in total)]
    if(find(~isfinite(ecgsig))) %check if the data is finite
        beats_inflation{i} = [];
%         winStart_infl(i) = 0;
%         winEnd_infl(i) = 0;
        disp(['WARNING: Patient ', num2str(i), ' ignored due to non-finite elements']);
        continue;
    end
    
    %Check the minimum number of samples for inflation patients - 55000 (55 seconds)
    N_min = min(inflationStartEnd_n(2,:) - inflationStartEnd_n(1,:));
    trim = 1000; %ignore the first and last second
    ecgsigInflation = ecgsig(inflationStartEnd_n(2,i)-N_min+trim : ...
                            inflationStartEnd_n(2,i)-trim, :);
    
    %% Average complexes (plot graphs for one patient)
    if((i == PATIENT_DEBUG) && DEBUG)
%         [beats_inflation{i},winStart_infl(i),winEnd_infl(i)] = averageBeat(ecgsigInflation, MIN_GAP, MIN_BEAT_LENGTH,1);
        beats_inflation{i} = averageBeat(ecgsigInflation,1);

    else
%         [beats_inflation{i},winStart_infl(i),winEnd_infl(i)] = averageBeat(ecgsigInflation, MIN_GAP, MIN_BEAT_LENGTH,0);
        beats_inflation{i} = averageBeat(ecgsigInflation,0);
    end
    disp(['Patient ', num2str(i), ' completed'])
end
clear N_min ecgSigAll i trim ecg ecgsig ecgsigInflation N_Min ...
    PATIENT_DEBUG 
disp('**** Inflation complexes complete ****');
disp('');

%% Package average beats with annotations
% Find empty cells
blanks_base = find(cellfun(@isempty,beats_baseline));
blanks_infl = find(cellfun(@isempty,beats_inflation));
%delete empty cells and filenames
beats_baseline = beats_baseline(~cellfun('isempty',beats_baseline));
beats_inflation = beats_inflation(~cellfun('isempty',beats_inflation));
filesBaseline(blanks_base) = [];
filesInflation(blanks_infl) = [];
% %delete empty window start/ends
% winStart_base(blanks_base) = [];
% winEnd_base(blanks_base) = [];
% winStart_infl(blanks_infl) = [];
% winEnd_infl(blanks_infl) = [];

% Get annotation row for each file (baseline/inflation)
for i = 1:length(filesBaseline)
    strB{i}= filesBaseline{i}(7:end-4);  %find filename
    %find row in ann containing filename and store in ann
    ann_baseline(i,:) = ann(strcmp(ann.filename,strB{i}),:);
end
ann_inflation = array2table([]);
h_old=0; h_new=0;
for i = 1:length(filesInflation)    %for all inflation beats
    strI{i}= filesInflation{i}(7:end-4);  %find filename
    %find row in ann containing filename and store in ann
%     ann_inflation(i,:) = ann(strcmp(ann.filename,strI{i}),:);
    if(h_new-h_old > 1) %check if two strings of the same name have been found
        h_old = h_old + 1;  
        continue;   %avoid duplication of rows by skipping current row
    end
    h_old = height(ann_inflation);  %check the current height of the table
    % concatenate the rows of inflation patients in ann to ann_inflation
    ann_inflation = [ann_inflation; ann(strcmp(ann.filename,strI{i}),:)];  
    h_new = height(ann_inflation);  %check the new height of the table to avoid dupes
end
clear blanks_base blanks_infl strI strB h_new h_old

%% Split inflation files into LAD,LCX,RCA
%find and extract all LAD rows
ann_inflation.artery = string(ann_inflation.artery); %convert categorical to string
ann_LAD = ann_inflation(contains(ann_inflation.artery,"LAD"),:);
ann_LCX = ann_inflation(contains(ann_inflation.artery,"circ"),:);
ann_RCA = ann_inflation(contains(ann_inflation.artery,"RCA"),:);

%% Split average inflation beats into LAD,LCX,RCA
beats_LAD = beats_inflation(contains(ann_inflation.artery,"LAD"));
beats_LCX = beats_inflation(contains(ann_inflation.artery,"circ"));
beats_RCA = beats_inflation(contains(ann_inflation.artery,"RCA"));

%% Run STEMI detector - baselines
disp("Starting STEMI detection");
%get input arguments
jPoint_n_base = ones(1,length(filesBaseline)) * 470; %set all jpoints to 470 temporarily
ann_baseline.gender = string(ann_baseline.gender); %convert categorical to string
maleFemale_base = strcmpi(ann_baseline.gender,"F"); %find females ;) (1)
age_base = ann_baseline.age; %extract the baseline ages

result_base = zeros(1,length(filesBaseline));   %preallocation
location_base = zeros(1,length(filesBaseline)); %preallocation
for i=1:length(filesBaseline)
    [result_base(i), location_base(i)] = detectSTEMI(beats_baseline(i), jPoint_n_base(i), maleFemale_base(i), age_base(i));
end

%% Run STEMI detector - inflation
jPoint_n_inf = ones(1,length(filesInflation)) * 470; %set all jpoints to 470 temporarily
ann_inflation.gender = string(ann_inflation.gender); %convert categorical to string
maleFemale_inf = strcmpi(ann_inflation.gender,"F"); %find females ;) (1)
age_inf = ann_inflation.age; %extract the baseline ages

result_inf = zeros(1,length(filesInflation));   %preallocation
location_inf = zeros(1,length(filesInflation)); %preallocation
for i=1:length(filesInflation)
    [result_inf(i), location_inf(i)] = detectSTEMI(beats_inflation(i), jPoint_n_inf(i), maleFemale_inf(i), age_inf(i));
end
disp("STEMI detection complete");

%% Plot a single patient
if(DEBUG)
    i = 3;
    ecgsig = importdata(filesInflation(i));  %import all leads for annotation (V1-6; I-III; 9 in total)]
    if(find(~isfinite(ecgsig))) %check if the data is finite
        disp(['WARNING: Patient ', num2str(i), ' should be ignored due to non-finite elements']);
    end
    %Check the minimum number of samples for inflation patients - 55000 (55 seconds)
    N_min = min(inflationStartEnd_n(2,:) - inflationStartEnd_n(1,:));
    trim = 1000; %ignore the first and last second
    ecgsigInflation = ecgsig(inflationStartEnd_n(2,i)-N_min+trim : ...
                            inflationStartEnd_n(2,i)-trim, :);
    [~] = averageBeat(ecgsigInflation,1);
end

%Improvements for future:
% - Dynamic N_DELAY for varying resting heart rates
% - Window works backwards from the deflation time (if applicable)
%    instead of a small 53 second window at the end
% - Get average complex for pre-inflation alongside PBI like horacek data