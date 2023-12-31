%% Plot BSPM-based leads against each other in the time domain
function compareLeads(lead1,lead2,fs)
    N = length(lead1(1,:));
    t = [0:N-1]*(1/fs);
%     set(gcf,'position', [100, 400, 300, 300]);   %plot is 300x300 wide
    set(gcf,'position', [100, 100, 300, 600]);   %plot is 300x300 wide
    % hax=axes; 
    plot(t,lead1, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':');
    hold on;
    plot(t,lead2, 'LineWidth', 1.5, 'Color', 'k');
    title('Measured Lead vs Augmented Lead');
    xlabel('Time (s)');
    ylabel('Amplitude ($\mu$V)');
%     SP=t(result{1,pno_result}(1,5)); %your point goes here 
%     xline(SP, 'Color',[0 0 0], 'LineWidth', 1, 'LineStyle', '--');
    legend('Measured','Augmented','Location','southeast');
%     xlim([t(40) t(350)]); %control how many points to plot (narrow x axis)
    hold off; grid on;
end