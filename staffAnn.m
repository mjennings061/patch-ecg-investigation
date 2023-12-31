%% Manual annotation of STAFF data in SSL format
clc;clear;close all;

%% Import median beat STAFF data
load('beats_baseline.mat'); %import baseline median complexes
load('beats_inflation.mat');
beats = [beats_baseline beats_inflation];
clear beats_baseline beats_inflation;

%% Get coefficients for 8-lead to SSL
load('coeffs.mat'); %8-lead to SSL coeffs

%% Generate SSLS
beats8 = {}; %8 lead beats
ssl_ste = {}; %short spaced STE sensitive lead
order = [7 8 1 2 3 4 5 6];
for i=1:length(beats)
    %format STAFF data to 8-lead format
    for j=1:8 %each lead of the new 8-lead format
        beats8{i}(:,j) = beats{i}(:,order(j));
    end
    
    %perform regression calculation
    ssl_ste{i} = c_steLead * beats8{i}';    %new lead (transformed)
    ssl_orth{i} = c_orthLead * beats8{i}';    %new lead (transformed)
end

%% Plot and annotate
figure('Position', [800 100 700 400]);
%plot complex
%ask for start of P
%ask for end of P
%ask for start of QRS
%ask for end of QRS
%ask for start of T
%ask for end of T

%% Repackage files like patchMI.m (Nx6 data; Nx1 annotations)

%package and save as cell
in = input('Save SSLs? (y/n)','s');
if(in == 'y')
    save('ssls', 'ssl_orth', 'ssl_ste');
    disp('staffAnn: Short spaced leads saved as ssls.mat');
end
%% Clear variables
clear i in order j;