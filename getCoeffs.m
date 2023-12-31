%% Generate transform coefficients from 12-lead to a selected lead
%Inputs:
% - data - 352 node BSPM data
% - node_p - Selected +ve node number on dalhousie torso pdf file
% - node_n - Selected -ve node number on dalhousie torso pdf file
%Output:
% - coeffs - A set of 1x8 transform coefficients for the chosen lead
function coeffs = getCoeffs(data, node_p, node_n)
tic;
    for i = 1:length(data)  %for all patients
        %% Extract real lead data
         %One matrix for all leads wrt time
        %         time (t) -->-->
        % lead I [n n-1 n-2 ... n-N]
        % lead II[n n-1 n-2 ... n-N]
        % ...
        % lead V6[n n-1 n-2 ... n-N]
        N = length(data{i}(1,:));
        realLead = zeros(8,N);      %row = lead I to V6. col = samples(1:N)
        targetLead = zeros(1,N);   %row = chosen lead. col = samples(1:N)
    %     RA = (data{pNum}(63,:) + data{pNum}(104,:))/2; %ML config
    %     LA = (data{pNum}(53,:) + data{pNum}(93,:))/2; %ML config
    %     LL = (3*data{pNum}(346,:) + 2*data{pNum}(347,:))/5; %ML config
    %     realLead(1,:) = LA-RA; %lead I
    %     realLead(2,:) = LL-RA; %lead II
        realLead(1,:) = data{i}(2,:); %I
        realLead(2,:) = data{i}(3,:); %II
        realLead(3,:) = data{i}(172,:); %V1
        realLead(4,:) = data{i}(174,:); %V2
        realLead(5,:) = (data{i}(195,:) + data{i}(196,:))/2; %V3
        realLead(6,:) = data{i}(219,:); %V4
        realLead(7,:) = (data{i}(220,:) + 2*data{i}(221,:))/3; %V5
        realLead(8,:) = data{i}(222,:); %V6
        targetLead(1,:) = data{i}(node_p+3,:) - data{i}(node_n+3,:);

        %% Calculate regression coefficient
        % B = Y/X in Kartheeban Nagenthiraja's work
        % where B = corr coeffs, Y = V7-V12 and X = lead I-V6
       B(i,:) = targetLead(1,:)/realLead;    %V7 coefficients for each patient
    end
    coeffs = zeros(1,8);
    for i=1:8
        coeffs(1,i) = median(B(:,i));
    end
t = toc;
disp(['getCoeffs: ', num2str(t), ' seconds']);
end