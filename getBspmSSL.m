%% Transform 12-lead ECG from BSPM data using coeffs. Return root mean square error and correlation coefficient
%Inputs:
% - bspmData - 352 node BSPM data for a single patient. In matrix format, not cell format
% - coeffs - A set of 1x8 transform coefficients for the chosen lead
% - node_p - Selected +ve node number on dalhousie torso pdf file
% - node_n - Selected -ve node number on dalhousie torso pdf file
%Output:
% - leadOut - Time domain signal of the transformed lead
% - rmse - Root mean square error between chosen lead and the actual 
function [leadOut, rmse, cc] = getBspmSSL(bspmData, coeffs, node_p, node_n)
    N = length(bspmData(1,:));
    lead = zeros(8, N); %row=lead I to V6. col = samples(1:N)
%     RA = (data{pNum}(63,:) + data{pNum}(104,:))/2;
%     LA = (data{pNum}(53,:) + data{pNum}(93,:))/2;
%     LL = (3*data{pNum}(346,:) + 2*data{pNum}(347,:))/5;
%     lead(1,:) = LA-RA; %lead I
%     lead(2,:) = LL-RA; %lead II
    lead(1,:) = bspmData(2,:); %I
    lead(2,:) = bspmData(3,:); %II
    lead(3,:) = bspmData(172,:); %V1
    lead(4,:) = bspmData(174,:); %V2
    lead(5,:) = (bspmData(195,:) + bspmData(196,:))/2; %V3
    lead(6,:) = bspmData(219,:); %V4
    lead(7,:) = (bspmData(220,:) + 2*bspmData(221,:))/3; %V5
    lead(8,:) = bspmData(222,:); %V6
    leadOutReal = bspmData(node_p+3,:) - bspmData(node_n+3,:);  %actual lead
    
    leadOut = coeffs * lead;    %new lead (transformed)
    
    %% Calculate RMSE and Correlation Coefficient
    rmse = sqrt(mean((leadOutReal(1,:) - leadOut(1,:)).^2));  % Root Mean Squared Error in uV
    c = corrcoef(leadOutReal(1,:),leadOut(1,:)); %correlation coefficient (0 to 1)
    cc = c(1,2); %extract cc from 2x2 matrix (off-diagonals e.g 1,2;2,1)
end