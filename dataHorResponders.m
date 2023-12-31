%% Filter horacek data for responders only
%Input:
% fileNames - Patient information in a txt file from the horacek dataset
% vessel - Vessel to filter e.g. 'LAD','LCX','RCA'
%Output:
% patients - positive responding patient names
% patientIdx - the index in BalloonBSPMdata of positive responders
function [patients, patientIdx] = dataHorResponders(fileNames,vessel)
tic;
    patients = {};
    patientIdx = [];
    j=1;
    if(nargin == 1)
        %% All responding patients 
        for i=1:length(fileNames) 
            if (contains(fileNames(i),'_Y_')) %Only look at responders
                if(contains(fileNames(i), '_B_'))
                    if(contains(fileNames(i+1),'_P_'))
                        patientIdx(j,1) = i; %patient number is stored for future reference
                        patients(j,1) = fileNames(i);
                        j=j+1;
                        patientIdx(j,1) = i+1;
                        patients(j,1) = fileNames(i+1);
                        j=j+1;
                        i=i+2;
                    end
                end
            end
        end
    elseif(nargin == 2)
        %% Responding patients in one vessel
        for i=1:length(fileNames) 
            if (contains(fileNames(i),'_Y_') && contains(fileNames(i),vessel)) %Only look at responders
                if(contains(fileNames(i), '_B_'))
                    if(contains(fileNames(i+1),'_P_'))
                        patientIdx(j,1) = i; %patient number is stored for future reference
                        patients(j,1) = fileNames(i);
                        j=j+1;
                        patientIdx(j,1) = i+1;
                        patients(j,1) = fileNames(i+1);
                        j=j+1;
                        i=i+2;
                    end
                end
            end
        end
    end
t = toc;
disp(['dataHorResponders: ', num2str(t), ' seconds']);
end