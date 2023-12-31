%%Take already ranked data and filter out those above MAX_LEAD_DISTANCE in mm
%Inputs:
% - points - 3D torso data as a .pts file
% - face - 3D facet data as a .fac file
% - rankedData - ranked leads (3 cols: col 1 = electrode+; col 2 = electrode-; col 3 = rank)
% - MAX_LEAD_DISTANCE - the furthest apart the short spaced/vector lead can be in mm
%Output:
% - rankedDataSVL - All leads ranked below XXX mm in the same format as rankedData

function rankedDataSVL = filterLeadLength(points, face, rankedData, MAX_LEAD_DISTANCE)
tic;
    x = points(:,1);    %split off x,y,z coordinates 
    y = points(:,2);
    z = points(:,3);
    tr = triangulation(face,x(:),y(:),z(:));    %Create a surface plot with triangulation
    rankedDataSVLdupe = []; %blank array. This will hold the short leads including duplicate leads
    for i=1:length(rankedData(:,1))    %for all possible leads
        node_p = rankedData(i,2);   %node+ is column 2
        node_n = rankedData(i,1);   %node- os column 1
        if((node_p<=352) && (node_n<=352))
            D = pdist([points(node_p,:); points(node_n,:)]);    %calculate node distance
            if (D <= MAX_LEAD_DISTANCE)     %only allow those below the max distance
                rankedDataSVLdupe = [rankedDataSVLdupe; rankedData(i,:)];   %concatenate matrix with the next lead under XXX mm
            end
        end
    end
    
    %% filter out duplicate leads
    remove = zeros(length(rankedDataSVLdupe(:,1)),1);
    node_currP = 1;
    node_currN = 1;
    for i = 1:length(rankedDataSVLdupe)
        node_P = rankedDataSVLdupe(i,1);
        node_N = rankedDataSVLdupe(i,2);
        for j = 1:length(rankedDataSVLdupe)
            if(rankedDataSVLdupe(j,1) == node_N && rankedDataSVLdupe(j,2) == node_P) %if the inverse pattern has been seen before
                if(remove(i) == 0)
                    remove(j) = 1; %remove lead j
                end
            end
        end
    end
    rankedDataSVL = rankedDataSVLdupe(~remove,:); %duplicate leads removed by a logical index of 1
t = toc;
disp(['filterLeadLength: ', num2str(t), ' seconds']);
end