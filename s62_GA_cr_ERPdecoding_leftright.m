clc;clear  
close all

%% parameters
pp2do               = [1:23];

nsmooth             = 50;

plotSingleSubs      = 0;
colormap2use        = fliplr(brewermap(100, 'RdBu'));

%% load
s = 0;
for pp = pp2do;
    s = s+1;
    param = getSubjParam(pp);
    
    disp(['getting data from ' param.subjName]);
    
    load([param.path, 'saved_data/decoding_leftright', '__' param.subjName], 'dec');
    
    if nsmooth > 0
        factorOf1000hz = (1000/(256/2)); % after epoching, 1000/256. After decoding another times /2 (currently). 
        for x1 = 1:size(dec.correlation,1);
            dec.correlation(x1,:) = smoothdata(squeeze(dec.correlation(x1,:)),'movmean',round(nsmooth/factorOf1000hz));
            dec.eucdist(x1,:) = smoothdata(squeeze(dec.eucdist(x1,:)),'movmean',round(nsmooth/factorOf1000hz));
            dec.mahdist(x1,:) = smoothdata(squeeze(dec.mahdist(x1,:)),'movmean',round(nsmooth/factorOf1000hz));
       end
    end
    
    d1(s,:,:) = dec.correlation;
    d2(s,:,:) = dec.eucdist;
    d3(s,:,:) = dec.mahdist;

end

%% plot individuals?

if plotSingleSubs

figure('Name','CueValidity_v1_decoLvsR_singlesub','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 900])
for subj = 1:s    
    subplot(5,5,subj); hold on; title(pp2do(subj));
    plot(dec.time, squeeze(d3(subj,:,:)));
    plot(xlim, [50,50], '--k'); plot([0,0], ylim, '--k');
    xlim([dec.time(1), dec.time(end)]);
    ylim([40 65]);
end
legend(dec.labels);
saveas(gcf, 'CueValidity_v1_decoLvsR_singlesub_s23.jpg');

end

%% put back into structure
dec.correlation     = squeeze(mean(d1));
dec.eucdist         = squeeze(mean(d2));
dec.mahdist         = squeeze(mean(d3));

%% plot ---------------------------------------------------------------------------------------


%% ERP timecourses in selected channels
figure('Name','CueValidity_v1_decoLvsR_valid4con_overlay','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 400])
subplot(1,3,1); hold on; title('correlation');
[m_plot1] = frevede_errorbarplot(dec.time, d1(:,1,:), [0, 0, 1], 'se');
[m_plot2] = frevede_errorbarplot(dec.time, d1(:,2,:), [1, 0, 0], 'se');
[m_plot3] = frevede_errorbarplot(dec.time, d1(:,3,:), [1, 0, 1], 'se');
[m_plot4] = frevede_errorbarplot(dec.time, d1(:,4,:), [0, 1, 1], 'se');
plot(xlim, [50,50], '--k'); plot([0,0], ylim, '--k');
legend([m_plot1, m_plot2, m_plot3, m_plot4] , dec.labels);
xlim([dec.time(1), dec.time(end)]);
ylim([40 60]);

subplot(1,3,2); hold on; title('euclid. dist');
[m_plot1] = frevede_errorbarplot(dec.time, d2(:,1,:), [0, 0, 1], 'se');
[m_plot2] = frevede_errorbarplot(dec.time, d2(:,2,:), [1, 0, 0], 'se');
[m_plot3] = frevede_errorbarplot(dec.time, d2(:,3,:), [1, 0, 1], 'se');
[m_plot4] = frevede_errorbarplot(dec.time, d2(:,4,:), [0, 1, 1], 'se');
plot(xlim, [50,50], '--k'); plot([0,0], ylim, '--k');
legend([m_plot1, m_plot2, m_plot3, m_plot4] , dec.labels);
xlim([dec.time(1), dec.time(end)]);
ylim([40 60]);

subplot(1,3,3); hold on; title('mahalanobis dist');
[m_plot1] = frevede_errorbarplot(dec.time, d3(:,1,:), [0, 0, 1], 'se');
[m_plot2] = frevede_errorbarplot(dec.time, d3(:,2,:), [1, 0, 0], 'se');
[m_plot3] = frevede_errorbarplot(dec.time, d3(:,3,:), [1, 0, 1], 'se');
[m_plot4] = frevede_errorbarplot(dec.time, d3(:,4,:), [0, 1, 1], 'se');
plot(xlim, [50,50], '--k'); plot([0,0], ylim, '--k');
legend([m_plot1, m_plot2, m_plot3, m_plot4] , dec.labels);
xlim([dec.time(1), dec.time(end)]);
ylim([40 60]);

saveas(gcf, 'CueValidity_v1_decoLvsR_valid4con_overlay_s23.jpg');


