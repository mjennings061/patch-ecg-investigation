%% Find the patient index with the highest STE
% Inputs:
%   - BalloonBSPMdata - 352 node ecg data
%   - patientNo - indexes of patients from BalloonBSPMdata
%   - jDelay - Number of samples after the j-point annotation to get ST segment
%   - node_p - Positive electrode
%   - node_n - Negative electrode 
% Output:
%   - idxOut - 1xN matrix of indexes. Highest positive ST difference across SSL first
function idxOut = steMaxIdx(BalloonBSPMdata, patientNo, jDelay, node_p, node_n)
tic;
    stValue = [];
    %% loop through patients for the highest STE
    N = length(patientNo);
    for i = 1:2:N   %loop through each patient and find STE on chosen lead
        % find first j point
        for j = 1:length(BalloonBSPMdata{i}(1,:))
            if(BalloonBSPMdata{i}(1,j) == 3)
                jPoint = j+jDelay;
                break;
            end
        end
        % find second j point
        for j = 1:length(BalloonBSPMdata{i+1}(1,:))
            if(BalloonBSPMdata{i+1}(1,j) == 3)
                jPoint_2 = j+jDelay;
                break;
            end
        end
        stValue(i) = BalloonBSPMdata{i}(node_p+3, jPoint)-BalloonBSPMdata{i}(node_n+3, jPoint);
        stValue(i+1) = BalloonBSPMdata{i+1}(node_p+3, jPoint_2)-BalloonBSPMdata{i+1}(node_n+3, jPoint_2);
        ste(i) = stValue(i+1)-stValue(i);
    end
%     [maxSTE,idxOut] = max(ste);
    [maxSTE,idxOut] = sort(ste,'descend');
t = toc;
disp(['steMaxIdx: ', num2str(t), ' seconds']);
end