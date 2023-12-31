%% Calculate all possible leads and rank them based on ST elevation
function rankedData = steRank(BalloonBSPMdata, patientNo)
tic;
    %% Calculate STE for all possible patients and all lead combinations
    for h = 1:2:length(patientNo)-1
        PATIENT_NO = patientNo(h,1);
        c1 = BalloonBSPMdata{1,PATIENT_NO}(:,:); %Baseline data
        c2 = BalloonBSPMdata{1,PATIENT_NO+1}(:,:);   %Peak balloon data
        waveType = c1(1,:); %waveType is the annotation data for P QRS and ST
        waveType_2 = c2(1,:); %_2 is the annotation data for the PBI waveform
        for i=1:3
            c1(1,:) = [];   %remove the annotation data and LA RA electrodes
            c2(1,:) = [];
        end
        noNodes = length(c1(:,1));  %determine the number of total nodes (352)

        fs = 500; %sampling freq
        jDelay = 40e-3/(1/fs);  %a 40ms delay based after the annotated jpoint

        %% Find J-point on baseline and PBI samples
        %baseline
        for i=1:length(waveType)    %for all baseline annotations
            if(waveType(i) == 3)        %find the first j point annotation
                jPoint = i+jDelay;          %add the 40ms delay
                break
            end
        end

        %PBI
        for i=1:length(waveType_2)  %for all peak balloon inflation annotations
            if(waveType_2(i) == 3)
                jPoint_2 = i+jDelay;
                break
            end
        end

        %% Calculate all leads
        % Lead 1-1, lead 1-2
        % Baseline data only
        for i = 1:noNodes
            leadMatrix_b{i} = zeros(noNodes, length(c1(1,:))); %preallocate a 352xSampleLength matrix for all nodes
        end
        for i = 1:noNodes
            for j = 1:noNodes
                leadMatrix_b{i}(j,:) = c1(i,:) - c1(j,:);   %take the first node away from the second node to create a new lead
            end
        end

        % Lead 1-1, lead 1-2
        % Peak-balloon inflation data only
        for i = 1:noNodes
            leadMatrix_p{i} = zeros(noNodes, length(c2(1,:)));
        end
        for i = 1:noNodes
            for j = 1:noNodes
                leadMatrix_p{i}(j,:) = c2(i,:) - c2(j,:); %peak balloon inflation leads
            end
        end

        % leadMatrix compares baseline (cell row 1) to PBI (cell row 2)
        % e.g. leadMatrix{1,2}(3,:) = {baselineData, node2}(node2-node3, :allSamples)
        leadMatrix = [leadMatrix_b; leadMatrix_p];

        %% Calculate ST-elevation at j-point for all possible leads
        %for each lead in leadMatrix_b
        for i=1:noNodes
            leadMatrix{3,i} = zeros(noNodes,7); %make array for st-elevations
        end

        %Find the ST elevation in all possible leads
        for i=1:noNodes %for all nodes (1-352)
            for j=1:noNodes %for all leads 352*(1-352)
                leadMatrix{3,i}(j,1) = i; %positive node number
                leadMatrix{3,i}(j,2) = j; %negative node number
                for k=1:length(waveType) %for all samples
                    if(waveType(k) == 3)
                        jPoint = k+jDelay; %Delay by jDelay samples (20 or 40ms)
                        leadMatrix{3,i}(j,3:4) = [jPoint leadMatrix{1,i}(j,jPoint)]; %Log the j point and amplitude at baseline
                        break;
                    end
                end
                for k=1:length(waveType_2)
                    if(waveType_2(k) == 3)
                        jPoint_2 = k+jDelay;
                        leadMatrix{3,i}(j,5:6) = [jPoint_2 leadMatrix{2,i}(j,jPoint_2)]; %Log the j point and amplitude at PBI
                        break;
                    end
                end
                leadMatrix{3,i}(j,7) = leadMatrix{3,i}(j,6) - leadMatrix{3,i}(j,4); %log the ST elevation between PBI/baseline
            end
        end

        leadMatrix{4,1} = [];
        for i=1:noNodes
            leadMatrix{4,1} = [leadMatrix{4,1}; leadMatrix{3,i}];   %Put all data into one matrix
        end

        %% Sort and rank ST elevation at the j-point from highest to lowest
        result{1,h} = leadMatrix{4,1}; %Store current calculations into a large array for each patient
        result{2,h} = sortrows(result{1,h}, -7, 'ComparisonMethod', 'abs'); %Sort st elevation from highest to lowest
        result{1,h}(:,8) = 0;

        for rank=1:length(result{2,h})
            index = (result{2,h}(rank,1)-1)*noNodes+result{2,h}(rank,2);
            result{1,h}(index,8) = rank;    %Put the rank into the result matrix
        end

        if(h==1)
            result{3,1} = [result{1,1}(:,1:2)]; %put all ranked data into one array
            result{3,1}(:,3) = 0;
        end
        result{3,1}(:,3) = result{3,1}(:,3) + result{1,h}(:,8);  
    end
    
    rankedData = sortrows(result{3,1}, 3); %sort all combined nodes by biggest difference
    
t = toc;
disp(['steRank: ', num2str(t), ' seconds']);
end