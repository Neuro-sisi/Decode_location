%% Decoding Cueside from Alpha
clc;clear
close all

cd '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/EEGdata_Ana_FvE/'

%% parameters
for pp          = 1:23;

    plotResults     = 1;
    colormap2use    = fliplr(brewermap(100, 'RdBu'));

    %% load data
    param = getSubjParam(pp);

    load([param.path, '/processed_data/' 'epoched_data_eeg' '__', param.subjName], 'data');

    %% keep only channels of interest
    cfg         = [];
    cfg.channel = {'EEG'};
    data        = ft_preprocessing(cfg, data);

    %% remove bad ICA components
    load([param.path, '/saved_data/' 'ICAcomponents', '__' param.subjName], 'comp2rem','ica');
    cfg             = [];
    cfg.component   = comp2rem;
    data            = ft_rejectcomponent(cfg, ica, data);

    %% remove bad trials
    load([param.path, '/processed_data/' 'usableTrials', '__' param.subjName], 'trl2keep');

    cfg         = [];
    cfg.trials  = trl2keep;
    data        = ft_selectdata(cfg, data);

    %% baseline correct
    cfg                 = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-.25 0]; % 250 ms pre-cue baseline
    data                = ft_preprocessing(cfg, data);

    % lpfilter
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 30;
    data         = ft_preprocessing(cfg, data);

    %% restrict range to speed up
    cfg         = [];
    cfg.latency = [-0.25 1.5];
    data        = ft_selectdata(cfg, data);

    %% resample further to speed up
    resamplefactor = 2;
    data.fsample = 256 ./ resamplefactor;
    for trl = 1:length(data.trial)
        data.time{trl} = data.time{trl}(1:resamplefactor:end);
        data.trial{trl} = data.trial{trl}(:,1:resamplefactor:end);
    end

    %% loop over 4 conditions
    for condition = 1:4 % 1 = high100, 2 = med80, 3 = low60, 4 = imper100

        % start from clear variables
        corr_correctclassification = [];
        dist_correctclassification = [];
        mahdist_correctclassification = [];

        %% get data into usable matrix of trl x elec x time
        cfg             = [];
        cfg.trials      = data.trialinfo(:,2)==condition;
        cfg.keeptrials  = 'yes';
        tl              = ft_timelockanalysis(cfg, data); % data in tl.trial matrix
        tl.time         = tl.time*1000; % s to ms

        left     = ismember(tl.trialinfo(:,1), [41]);
        right    = ismember(tl.trialinfo(:,1), [42]);

        % alpha calculation--hilbert transform
        freqs = [8 12];
        % low-pass filtering
        filtData = nan(size(tl.trial)); % trl*elec*time
        % get alpha power
        for c = 1:size(tl.trial,2) %elec
            filtData(:,c,:) =abs(hilbert(eegfilt(squeeze(tl.trial(:,c,:)),data.fsample,freqs(1,1),freqs(1,2))')').^2;   %Instantaneous power
        end
        tl.trial = filtData; clear filtData


        %% get decoding per condition
        %% loop over trials and timepoints to get similarity to matching and non-matching classes (i.e. the actual "decoding"), and calculate decoding accuracy
        % using a "leave-one-out" approach (i.e., iterative "testing" each trial against all remaining trials)
        % logic: if pattern distinguishes between the two classes, the pattern in a given trial (at a given timepoint) should be MORE SIMILAR to other trials from the MATCHING (same) vs. NON-MATCHING (different) class.
        % => we therefore just need to look at pattern similarities (below, we quantify similarity in three different ways).

        d = tl.trial; % call data matrix just d, for simplifying script below -- d = trials x electrodes x timepoints

        alltrials = 1:size(d,1); % 1:ntrials
        classes = left; % 1 for this condition; 0 for the alternative condition (2 class decoding)

        for trl = alltrials % loop over all trials
            disp(['decoding trial ', num2str(trl), ' out of ', num2str(length(alltrials)), ' -- condition ', num2str(condition), ' -- pp ' num2str(pp)]);

            % "test data" (the data from this trial)
            thistrial = trl;
            class_thistrial = classes(trl);
            testdata = squeeze(d(thistrial,:,:));

            % "training data" (the matching and non-matching data from all remaining trials)
            othertrials         = ~ismember(alltrials, thistrial); % mark which trials are NOT the current trial - to compare the current trial against (never include "test" data in "training set")
            other_matching      = othertrials' & (classes==class_thistrial); %  trials from same class to which current trial belongs (EXCEPT the one we are currently looping over - i.e. "leave-one-out")
            other_nonmatching   = othertrials' & (classes~=class_thistrial); % trials for class to which current trials does NOT belong
            traindata_match     = squeeze(mean(d(other_matching,:,:))); % average over all matching trials
            traindata_nonmatch  = squeeze(mean(d(other_nonmatching,:,:))); % average over all non-matching trials

            %     covar = cov(squeeze(d(trl,:,:))'); % get covariance across time, to be used below in 'mahalanobis'.
            %     covar = nearestSPD(covar); % ensure its always Positive Definite, even if data is reduced rank after interpolation and/or ICA

            % now loop over time points to run classification separately for each timepoint
            for time = 1:length(tl.time)

                %--% three different ways to "quantify" similarity to matching and non-matching classes of data (calculate and store separately for each trial and timepoint in the loop)

                % metric 1 - correlation of topographies to test "similarity" of pattern (multivariate pattern across electrodes) between current trial and all remaining matching vs. non-matching trials
                corr_with_match      = corr(testdata(:,time), traindata_match(:,time)); % how "similar" (here correlated) is current trial to mean of all other trials from SAME class?
                corr_with_nonmatch   = corr(testdata(:,time), traindata_nonmatch(:,time)); % how "similar" is current trial to mean of all other trials from DIFFERENT class?
                corr_correctclassification(trl,time) = corr_with_match > corr_with_nonmatch; % correct classification if LARGER correlation with MATCH trials... incorrect otherwise

                % metric 2 - euclidean distance in multi-dimentional electrode space
                dist_to_match       = pdist([testdata(:,time)'; traindata_match(:,time)'], 'euclidean');
                dist_to_nonmatch    = pdist([testdata(:,time)'; traindata_nonmatch(:,time)'], 'euclidean');
                dist_correctclassification(trl,time) = dist_to_match < dist_to_nonmatch; % correct classification if SMALLER distance to MATCH trials... incorrect otherwise

                % metric 3 - mahalanobis distance -- also taking covariance into account - i.e. distances normalised by the covariance structure in the data... should be more appropriate/sensitive
                covar = cov(squeeze(d(othertrials,:,time))); % covariance over all other trials, at this time point
                covar = nearestSPD(covar); % ensure its always Positive Definite, even if data is reduced rank after interpolation and/or ICA
                mahdist_to_match    = pdist([testdata(:,time)'; traindata_match(:,time)'], 'mahalanobis', covar);
                mahdist_to_nonmatch = pdist([testdata(:,time)'; traindata_nonmatch(:,time)'], 'mahalanobis', covar);
                mahdist_correctclassification(trl,time) = mahdist_to_match < mahdist_to_nonmatch; % correct classification if SMALLER distance to MATCH trials... incorrect otherwise

                %--% end loops over timepoints and trials
            end % end of loop over time points
        end % end of loop over trials

        %% get decoding score and save per condition
        dec_correlation(condition,:) = squeeze(mean(corr_correctclassification))*100;
        dec_eucdist(condition,:)     = squeeze(mean(dist_correctclassification))*100;
        dec_mahdist(condition,:)     = squeeze(mean(mahdist_correctclassification))*100;

    end % end condition loop

    %% put into structure to plot and save
    dec = [];
    dec.time = tl.time;
    dec.labels = {'imper100','high100','med80','low60'};
    dec.correlation = dec_correlation([4,1,2,3],:);
    dec.eucdist = dec_eucdist([4,1,2,3],:);
    dec.mahdist = dec_mahdist([4,1,2,3],:);

    %% plot
    if plotResults
        figure;
        subplot(1,3,1); hold on; plot(dec.time, dec.correlation); plot(xlim, [50 50], '--k');    legend(dec.labels); title('correlation');
        subplot(1,3,2); hold on; plot(dec.time, dec.eucdist); plot(xlim, [50 50], '--k');     ylim([25 75]);  legend(dec.labels); title('euclidean dist');
        subplot(1,3,3); hold on; plot(dec.time, dec.mahdist); plot(xlim, [50 50], '--k');     ylim([25 75]);  legend(dec.labels); title('mahalanobis dist');
    end

    drawnow;

    %% save
    save([param.path, 'saved_data/decoding_alpha_leftright', '__' param.subjName], 'dec');

    %% end pp loop
end