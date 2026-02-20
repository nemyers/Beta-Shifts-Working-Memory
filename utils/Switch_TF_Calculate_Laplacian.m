clearvars; close('all'); clc

workstation = 'nick_pc';
switch workstation
    case 'nick_pc'
        onedrivepath  = 'E:\Myers\OneDrive - The University of Nottingham';
        fieldtrippath = [onedrivepath '\Matlab\toolbox\fieldtrip\'];
        studypath     = [onedrivepath '\Projects\El-SwitchBeta'];        
end

cd(studypath)
addpath(fieldtrippath);
ft_defaults;
fs = filesep;
toolboxpath   = [onedrivepath fs 'Projects' fs 'toolbox'];
addpath(toolboxpath);
SetupToolboxPath(toolboxpath);
addpath([studypath fs 'scripts'])

resultspath = [studypath];

datapath = [studypath '/data/derivatives'];

sublist = [2 4 5 7 11 14 15 16 17 18 19 21 22 23 27 29 31 32 35 36 2:31];
explist = [ones(1,20)*1 ones(1,30)*2];
nsubs   = length(sublist);
subpairs = [05 10 11 12 13 14 15; 
            28 23 22 24 34 42 41];
%%
for isub = 1:nsubs %subject loop
    % load data
    experiment = explist(isub);    
    substrg  = sprintf('S%02d',sublist(isub));
    
    subpath  = sprintf('%s%sLaplace%sSwitch_EEG%d%s%s',datapath,fs,fs,experiment,fs,substrg);
    %eyetrackpath = sprintf('%seyetracking',studypath);
    filename = sprintf('%s%s%s_EEG_final_laplace.mat',subpath,fs,substrg);
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);
    ntrials = length(data.trial);
    
    try, data = rmfield(data,'nTrials'); catch; end
    time    = data.time{1};
    
    outpath = sprintf('%s%sTF_phase%sSwitch_EEG%d',datapath,fs,fs,experiment);
    if ~exist(outpath,'dir'); mkdir(outpath); end
    % setup TF
    %%  do TF analysis, or load if it already exists
    overwrite_TF = 1;
    fname     = sprintf('%s%s%s_HanningTF_Laplace.mat',outpath,fs,substrg);
    fname_phs = sprintf('%s%s%s_HanningTF_Laplace_phs.mat',outpath,fs,substrg);
    lock_to_target = false;
    if lock_to_target
        fname = sprintf('%s%s%s_HanningTF_Laplace_Probe.mat',outpath,fs,substrg);
    end
    if ~exist(fname,'file') | overwrite_TF
        % Preprocessing Settings
        chanlist = 1:61; %skip EOG channels (61 and 62)
        nchans   = length(chanlist);

        fsample  = data.fsample;    
        dtime    = fix([-1.000 +2.996]*fsample); %time window of interest, relative to stim1 onset
        if lock_to_target
            dtime    = fix([-2.000 +1.000]*fsample); %time window of interest, relative to cue onset
        end
        dstep    = 1;
        timelist = (dtime(1):dstep:dtime(2))/fsample;
        ntimes = length(timelist);

        % re-cut epochs around probe onset and possibly remove baseline
        %waitbar(0,hbar,sprintf('processing S%02d...',sublist(isub)));
        ntrials = length(data.trial);
        datanew = data;
        dataerp_elem = nan(ntrials,nchans,ntimes);
        for itrial = 1:ntrials
            ionset = find(data.time{itrial} <= 0,1,'last');
            if lock_to_target
                ionset = find(data.time{itrial} <= 1.800,1,'last');
            end
            itime  = ionset+(dtime(1):dstep:dtime(2));
            itime  = min(itime,size(data.trial{itrial},2));
            dataerp_elem(itrial,:,:) = data.trial{itrial}(chanlist,itime);
            datanew.trial{itrial} = squeeze(dataerp_elem(itrial,:,:));
            datanew.time{itrial}  = timelist;
        end
        data = datanew; 
        data.label = data.label(chanlist);
        clear datanew dataerp_elem

        % Add zero-padding
        %dont use this yet
        add_padding = 1; %use zero-padding to avoid TF edge artifacts?
        if add_padding
            pad = ceil([1.000 1.000]*fsample)/fsample; %1  'PO3' 'PO4' 'P3' 'P4'sec should be plenty for alpha (8 Hz = 125 ms, at 5 cycles per wavelet = 750 ms)
            nchans = size(data.trial{1},1);
            for itrial = 1:ntrials
                data.trial{itrial} = [zeros(nchans,pad(1)*fsample),data.trial{itrial},zeros(nchans,pad(2)*fsample)];
                data.time{itrial}  = linspace(data.time{itrial}(1)-pad(1),data.time{itrial}(end)+pad(2),size(data.trial{itrial},2));
            end
        end
        % Time-Frequency Analysis - Settings
        tfmeth = 'hanning'; % time-frequency analysis method (multitapering or hanning or wavelet)
        chanoi = {};

        tons = -0.500;
        toff = +1.600;
        if lock_to_target
            tons = -1.500;
            toff = +0.500;
        end

        timeoi = data.time{1}(data.time{1} >= tons & data.time{1} <= toff); % time samples of interest (s)
        dstep  = 10; %used to be 5!
        timeoi = timeoi(1:dstep:end);
        triallen = timeoi(end)-timeoi(1);%window length in sec

        minf = 1; maxf = 40;
        freqoi = minf:maxf;

        switch lower(tfmeth)        
            case 'multitapering' % time-frequency analysis via multitapering (k = 2*tw*fw-1)
                cfg            = [];
                cfg.method     = 'mtmconvol';
                cfg.output     = 'pow';
                cfg.keeptrials = 'yes';
                cfg.feedback   = 'none';
                cfg.channel    = chanoi;
                cfg.foi        = freqoi;
                cfg.toi        = timeoi;
                cfg.taper      = 'dpss';
                cfg.t_ftimwin  = 5./cfg.foi; % time window (s)
                cfg.tapsmofrq  = 0.25*cfg.foi; % frequency window (Hz)

            case 'hanning' % time-frequency analysis via Hanning tapers (fw = 1/tw)
                cfg            = [];
                cfg.method     = 'mtmconvol';
                cfg.output     = 'pow';
                cfg.keeptrials = 'yes';
                %cfg.feedback   = 'none';
                cfg.channel    = chanoi;
                cfg.foi        = freqoi;
                cfg.toi        = timeoi;
                cfg.taper      = 'hanning';
                cfg.pad        = 9;
                cfg.t_ftimwin  = 5./cfg.foi; % time window (s) - 5 cycles
                %cfg.trials     = 1:2;

            case 'wavelet' % time-frequency analysis via Morlet wavelets (fw = f/w and tw = 1/fw)
                cfg            = [];
                cfg.method     = 'wavelet';
                cfg.output     = 'pow'; % fourier to get complex values, pow to get power
                cfg.keeptrials = 'yes';
                cfg.feedback   = 'none';
                cfg.channel    = chanoi;
                cfg.foi        = freqoi;
                cfg.toi        = timeoi;
                cfg.width      = 5; % wavelet width

            otherwise
                error('invalid time-frequency analysis method!');            
        end
        save_power = false;
        if save_power
            datatfr = ft_freqanalysis(cfg,data); %clear data
            datatfr.powspctrm = single(datatfr.powspctrm);
            save(fname,'datatfr','cfg','-v7.3');
        end

        save_phase = true;
        if save_phase
            minf = 4; maxf = 30;
            freqoi = minf:maxf;
            switch lower(tfmeth)        
                case 'hanning' % time-frequency analysis via Hanning tapers (fw = 1/tw)
                    cfg            = [];
                    cfg.method     = 'mtmconvol';
                    cfg.output     = 'fourier';
                    cfg.keeptrials = 'yes';
                    %cfg.feedback   = 'none';
                    cfg.channel    = chanoi;
                    cfg.foi        = freqoi;
                    cfg.toi        = timeoi;
                    cfg.taper      = 'hanning';
                    cfg.pad        = 9;
                    cfg.t_ftimwin  = 5./cfg.foi; % time window (s) - 5 cycles
               %     cfg.trials     = 1:50;        
            end
            dataphs = ft_freqanalysis(cfg,data); %clear data
            dataphs.fourierspctrm = single(dataphs.fourierspctrm);
            save(fname_phs,'dataphs','cfg','-v7.3');
        end
    else
        fprintf('\nAlready done!');
    end
end