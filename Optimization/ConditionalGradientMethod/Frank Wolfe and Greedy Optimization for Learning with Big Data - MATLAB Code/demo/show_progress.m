function show_progress(progress,x,s,p,T,r)
    % visualize progress, x and the current s
    subplot(2,1,1);
    plot(progress,'LineWidth',2,'Color',[1,0,0]);
    xlim([0 T]);
    legend('objective value')
    xlabel('iteration t')
    
    subplot(2,1,2);
    bar(1:p,s,1,'FaceColor',[0,0.7,0.7],'EdgeColor',[0,0.7,0.7]);
    ylim([-r r]);
    hold on
    
    bar(1:p,x,0.7,'FaceColor',[0.2,0.2,0.5],'EdgeColor','none');
    hold off
    xlabel('coordinate')
    legend('s','x') % add legend
    drawnow
end