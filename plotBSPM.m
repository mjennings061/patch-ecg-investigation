%% Plot a BSPM of the top patient with SSL electrode locations marked
%Input:
% - BalloonBSPMdata - 352 node BSPM data
% - patientNo - the index number in BalloonBSPMdata representing each patient
% - points - .pts file of torso
% - face - .fac file of faces
% - rankedDataSVL - 3 col data with node_n (col1), node_p (col2) and rank (col3)
% - noLeads - the number of SSLs to plot

function plotBSPM(BalloonBSPMdata, patientNo, points, face, rankedDataSVL, jDelay, noLeads)
tic;
    x = points(:,1);    %split off x,y,z coordinates 
    y = points(:,2);
    z = points(:,3);
    tr = triangulation(face,x(:),y(:),z(:));    %Create a surface plot with triangulation
    
    %% Calculate median amplitude at j+40ms across entire bspm
    allSTpoints = zeros(352,length(patientNo)/2); %preallocation
    for i = 1:length(patientNo)/2 %for each peak balloon slice at the stPoint
        for j = 1:length(BalloonBSPMdata{patientNo(i)}(1,:))
            if(BalloonBSPMdata{patientNo(i)}(1,j) == 3)
                stPoint = j+jDelay; %st point of the current patient
                break;
            end
        end
        allSTpoints(:,i) = BalloonBSPMdata{patientNo(i)}(4:end,stPoint); %record the amplitude of all electrodes at the ST point
    end
    medSTpoint = median(allSTpoints,2);
    
    %% Plot BSPM
    figure();
    [~] = trisurf(tr,medSTpoint);
    set(gcf,'position', [700, 400, 300, 300]);   %plot is 300x300 wide
    title('BSPM at J+40ms (during inflation)'); 
    xlabel('x (mm)');
    ylabel('y (mm)');
    zlabel('z (mm)');
    set(gca,'linewidth',1.4);
    hold on;
    view(0,270); %rotate the graph
    grid off;
    
    %% Plot the SSL
    plot3(x(rankedDataSVL(1,1)), y(rankedDataSVL(1,1)), z(rankedDataSVL(1,1)),'.w','markersize',15);
    plot3(x(rankedDataSVL(1,2)), y(rankedDataSVL(1,2)), z(rankedDataSVL(1,2)),'.w','markersize',15);
    plot3([x(rankedDataSVL(1,1)) x(rankedDataSVL(1,2))], ... 
        [y(rankedDataSVL(1,1)) y(rankedDataSVL(1,2))], ...
        [z(rankedDataSVL(1,1)) z(rankedDataSVL(1,2))], 'w');
    
    %% Additional SSLs
    if(noLeads > 1)
        for i = 2:noLeads
            plot3(x(rankedDataSVL(i,1)), y(rankedDataSVL(i,1)), z(rankedDataSVL(i,1)),'x','markersize',10);
            plot3(x(rankedDataSVL(i,2)), y(rankedDataSVL(i,2)), z(rankedDataSVL(i,2)),'x','markersize',10);
            plot3([x(rankedDataSVL(i,1)) x(rankedDataSVL(i,2))], ... 
            [y(rankedDataSVL(i,1)) y(rankedDataSVL(i,2))], ...
            [z(rankedDataSVL(i,1)) z(rankedDataSVL(i,2))]);
        end
    end
    
    %% Plot the V-leads
    plot3(x(169), y(169), z(169)-1,'.k','markersize',10); %v1
    plot3(x(171), y(171), z(171)-1,'.k','markersize',10); %v2
    plot3(x(192)+12, y(192), z(192)-1,'.k','markersize',10); %v3
    plot3(x(216), y(216), z(216)-1,'.k','markersize',10); %v4
    plot3(x(218), y(218), z(218)-1,'.k','markersize',10); %v5
    plot3(x(219), y(219), z(219)-1,'.k','markersize',10); %v6

    %% Configure axis and colourbar
    set(findall(gcf,'-property','FontWeight'),'FontWeight','bold');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(findall(gcf,'-property','FontSize'),'FontSize',9);
    c = colorbar;   %colour legend at the side
    [~] = c.LineWidth;    %set linewidth of colour map to 1.5
    c.LineWidth = 1.4;
    c.FontSize = 9;
    c.TickLabelInterpreter = 'latex';
    c.Label.Interpreter = 'latex';
    c.Label.String = 'Amplitude ($\mu$V)';
    c.Label.Position = [-1,43.0001964569092,0];
    shading interp; %blend the lines to remove meshing
    colormap jet;
    c.FontWeight = 'bold';

    %% Narrow the graph margins for publications
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    grid on;
    hold off;
    clear ax c hax p;
    
t = toc;
disp(['plotBSPM: ', num2str(t), ' seconds']);
end