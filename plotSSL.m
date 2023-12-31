%% Plot the time-domain orthogonal SSL
function plotSSL(BalloonBSPMdata, idx, node_p, node_n, node_p2, node_n2, fs)
tic;
    %check which signal is longer
    if(length(BalloonBSPMdata{idx+1}(1,:)) > length(BalloonBSPMdata{idx}(1,:)))
        t = [0:length(BalloonBSPMdata{idx+1}(1,:))-1] * 1/fs;
    else t = [0:length(BalloonBSPMdata{idx}(1,:))-1] * 1/fs;
    end
    
    figure()
    subplot(2,1,1)
    set(gcf,'position', [400, 400, 300, 300]);   %plot is 300x300 wide
    % hax=axes; 
    y_base = BalloonBSPMdata{idx}(node_p+3,:) - BalloonBSPMdata{idx}(node_n+3,:);
    y_peak = BalloonBSPMdata{idx+1}(node_p+3,:) - BalloonBSPMdata{idx+1}(node_n+3,:);
    plot(t(1:length(y_base)),y_base, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':');
    hold on;
    plot(t(1:length(y_peak)),y_peak, 'LineWidth', 1.5, 'Color', 'k');
    title('ST-Elevation Specific SSL');
    xlabel('Time (s)');
    ylabel('Amplitude ($\mu$V)');
    % SP=t(result{1,pno_result}(1,5)); %your point goes here 
    % xline(SP, 'Color',[0 0 0], 'LineWidth', 1, 'LineStyle', '--');
    legend('Baseline','PBI','Location','southeast');
%     xlim([t(40) t(350)]); %control how many points to plot (narrow x axis)
    hold off; grid on; 
    
    subplot(2,1,2)
    set(gcf,'position', [400, 100, 300, 600]);   %plot is 300x300 wide
    % hax=axes; 
    y_base2 = BalloonBSPMdata{idx}(node_p2+3,:) - BalloonBSPMdata{idx}(node_n2+3,:);
    y_peak2 = BalloonBSPMdata{idx+1}(node_p2+3,:) - BalloonBSPMdata{idx+1}(node_n2+3,:);
    plot(t(1:length(y_base2)),y_base2, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':');
    hold on;
    plot(t(1:length(y_peak2)),y_peak2, 'LineWidth', 1.5, 'Color', 'k');
    title('Orthogonal SSL');
    xlabel('Time (s)');
    ylabel('Amplitude ($\mu$V)');
    % SP=t(result{1,pno_result}(1,5)); %your point goes here 
    % xline(SP, 'Color',[0 0 0], 'LineWidth', 1, 'LineStyle', '--'); %plot a vertical line)
    legend('Baseline','PBI');
%     xlim([t(40) t(350)]); %control how many points to plot (narrow x axis)
    grid on;

    % This narrows the margins
    % ax = gca;
    % outerpos = ax.OuterPosition;
    % ti = ax.TightInset; 
    % left = outerpos(1) + ti(1);
    % bottom = outerpos(2) + ti(2);
    % ax_width = outerpos(3) - ti(1) - ti(3);
    % ax_height = outerpos(4) - ti(2) - ti(4);
    % ax.Position = [left bottom ax_width ax_height];
    % grid on;
    hold off;
t = toc;
disp(['plotSSL: ', num2str(t), ' seconds']);
end