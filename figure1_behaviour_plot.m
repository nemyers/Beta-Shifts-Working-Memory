% This script loads pre-calculated behavioral data and generates Figure 1
% from the paper, showing:
%   - Psychometric curves with swap error modeling
%   - Precision, guess rate, and swap rate by trial type
%   - Regression coefficients for item distance and congruency
%   - Item distance effects on accuracy and RT
%
% KEY FINDING:
%   Switch costs in accuracy and RT scale with item distance.
%   Swap errors are more frequent on switch trials.
%
% INPUTS:
%   - results/Study12_TFPow_Beh_Decoding.mat
%     (contains behavioral data for all subjects)
%
% OUTPUTS:
%   - Figure 1 panels (saved as PNG, EPS, and PDF)
%   - Statistical results printed to console
%
% Author: Nicholas E. Myers
% Cleaned for GitHub: 2026
% For: Myers, Stokes, & Muhle-Karbe (2026) J. Neuroscience

clearvars;
close all;
clc;

%% 1. SETUP AND LOAD DATA

% Load configuration
paths = switchbeta_config();

% Subject information
sublist = [2 4 5 7 11 14 15 16 17 18 19 21 22 23 27 29 31 32 35 36 2:31];
explist = [ones(1,20)*1 ones(1,30)*2];
nsubs   = length(sublist);

% Subject pairs (participated in both experiments)
subpairs = [05 10 11 12 13 14 15; 
            28 23 22 24 34 42 41];

% Load pre-calculated behavioral data
fprintf('Loading behavioral data...\n');
data_file = fullfile(paths.results, 'Study12_TFPow_Beh_Decoding.mat');
if ~exist(data_file, 'file')
    error('Behavioral data file not found: %s\nRun figure1_behavior_calculate.m first.', data_file);
end
load(data_file);
nsub = length(tfdata);
%%
outdata = [];
dvdata = [];
repeffect = [];
b = '';
fprintf('\nAnalysing: ');
for isub = 1:nsub %subject loop
    m = sprintf('S%02d/%02d',isub,nsub);
    fprintf([b m]); b = repmat('\b',[1 length(m)]);
    beh        = tfdata(isub).behaviour;
    experiment = tfdata(isub).experiment;
    subind     = tfdata(isub).subind;
    explist(isub,1) = experiment;
    
    rsp = 2-beh.resp;
    rt  = beh.rt;
    ac  = beh.acc;
    dv  = round(beh.dv*1e5)/1e5;
    rep = beh.repnum;
    bt  = beh.blocktrial;
    bl  = beh.block;
    td  = circ_dist(beh.template_angles(:,1)*2*pi/180,beh.template_angles(:,2)*2*pi/180)/pi*90;
    td  = round(td*1e4)/1e4;
    td  = abs(td);
    
    udv = unique(dv(:,1));
    ndv = length(udv);
    
    utd = unique(td);
    ntd = length(utd);
    
    d = dv(:,1);
    id = abs(d)>45;
    d(id,1) = (90-abs(d(id,1))).*sign(d(id,1));
    
    dwrap = d;
    %ds = dv(:,1);
    %ids = ds>45;
    %ds(ids) = 90-ds(ids)
    d = abs(d);

    d2 = dv(:,2);
    id = abs(d2)>45;
    d2(id,1) = (90-abs(d2(id,1))).*sign(d2(id,1));
    dwrap(:,2) = d2;
    d2 = abs(d2);
    
    udv = unique(d);
    ndv = length(udv);
    dmin = -45; dmax = -dmin;

    
    congresp = diff(sign(beh.dv),[],2)==0;
    
    maxrep = 3;
    rep(rep>maxrep) = maxrep;
    
    blocktrialmin = 2;
    
    for irep = 1:maxrep
        itrl = bt>=blocktrialmin & rep==irep;
        repeffect(isub,irep,1) = mean(ac(itrl,1));
        itrl = itrl & ac==1;
        repeffect(isub,irep,2) = median(rt(itrl,1));        
    end
    %since performance is comparaable on repeat 1 and repeat 2+, average
    %all repeats
    %(using only repeat1 yields the same results)
    maxrep = 3;
    rep(rep>maxrep) = maxrep;
    for itd = 1:ntd
      for irep = 1:2
          itrl = bt>=blocktrialmin & rep==irep & td==utd(itd);

          %itrl = bt>=blocktrialmin & rep==irep & td==utd(itd) ;
          outdata(isub,itd,irep,1) = mean(ac(itrl),1);
          itrl = itrl & ac==1;
          outdata(isub,itd,irep,2) = median(rt(itrl),1);
      end
    end
    
    use_modelfit    =1;
    
    if(1)
    nd2 = length(unique(dwrap(:,1)));
    udv2 = unique(dwrap(:,1));

    
    udv2_abs = unique(abs(dwrap(:,1)));
    nd2_abs  = length(udv2_abs);
    for irep = 1:2
        for idv = 1:2
            for ind2 = 1:nd2
                itrl = bt>=blocktrialmin & rep==irep & dwrap(:,idv)==udv2(ind2); 
                if nnz(itrl)>0
                    misbnd_resp(isub,ind2,irep,idv)= mean(rsp(itrl));
                else
                    misbnd_resp(isub,ind2,irep,idv)= NaN;
                end
                itrl = itrl & ac==1;
                if nnz(itrl)>0
                    misbnd_rt(isub,ind2,irep,idv)= median(rt(itrl));
                else
                    misbnd_rt(isub,ind2,irep,idv)= NaN;
                end
            end
            for ind2_abs = 1:nd2_abs
                itrl = bt>=blocktrialmin & rep==irep & abs(dwrap(:,idv))==udv2_abs(ind2_abs); 
                if nnz(itrl)>0
                    misbnd_ac(isub,ind2_abs,irep,idv)= mean(ac(itrl));
                else
                    misbnd_ac(isub,ind2_abs,irep,idv)= NaN;
                end
            end
        end
    end
    end
           
    
    
     for idv = 1:ndv
        for irep = 1:2
            itrl = bt>=blocktrialmin & rep==irep & d(:,1)==udv(idv);
            dvdata(isub,idv,irep,1) = mean(ac(itrl,1));
        end
     end
     
    itrl = bt>=blocktrialmin;
    minrt = 0.1;
    maxrt = +Inf; maxrt = 2; 
    for itd = 1%:ntd
    for it = 1:2
        %itrl1 = rep == it & td == utd(itd) & d>dmin&d<dmax & itrl;
        itrl1 = rep == it & d>dmin&d<dmax & itrl;
        itrl2 = itrl1 & (ac == 1) & rt>minrt & rt<maxrt;
        tac(isub,it,itd) = nanmean(ac(itrl1));
        trt(isub,it,itd) = nanmedian(rt(itrl2));
        
        if use_modelfit
            design = [(d(itrl1))]; %absolute angular offset of the probe (in degrees)
            design = dwrap(itrl1,1);
            data   = [rsp(itrl1)]>0;
            %data   = [ac(itrl1)]>0;
            % create 2-parameter model with slope and lapse rate
            startp = [10.00 0.05 ]; %not sure if good starting values
            minp   = [0       0  ];
            maxp   = [+Inf    1 ];
            maxp   = [50      1 ];
            algopt     = optimset('fmincon');
            algopt     = optimset(algopt,'Algorithm','interior-point','Display','off'); %alternatives: 'active-set' or 'sqp' for fmincon
            [Pvec fval] = fmincon(@(P)costaccfun(P,design,data,false),startp,[],[],[],[],minp,maxp,[],algopt);
            mle_vals(isub,it,itd,:) = [Pvec fval];

            design = [d(itrl1) d2(itrl1)]; %absolute angular offset of the probe (in degrees)
            design = dwrap(itrl1,:);
            data   = [rsp(itrl1)]>0;
            %data   = [ac(itrl1)]>0;
            % create 2-parameter model with slope and lapse rate
            startp = [10.00 0.05 0.05]; %not sure if good starting values
            minp   = [0       0  0];
            maxp   = [+Inf    1 1];
            maxp   = [50      1 1];
            algopt     = optimset('fmincon');
            algopt     = optimset(algopt,'Algorithm','interior-point','Display','off'); %alternatives: 'active-set' or 'sqp' for fmincon
            [Pvec fval] = fmincon(@(P)costaccfun_misbind(P,design,data,false),startp,[],[],[],[],minp,maxp,[],algopt);
            mle_vals_misbind(isub,it,itd,:) = [Pvec fval];

            %predict responses
            respfun = @(P,x)((1-P(2)/2-P(3)).*normcdf(x(:,1),0,P(1)) + P(3).*normcdf(x(:,2),0,P(1)) + P(2)/2);
            yhat = respfun(Pvec,design);

            for idv = 1:2
                for ind2 = 1:nd2
                    misbnd_resp_hat(isub,ind2,it,idv)= mean(yhat(design(:,idv)==udv2(ind2)));
                end
            end
        end

        y = rsp(itrl1)>0;
        X = [dwrap];
        %X = [d ];
        X = X(itrl1,:);
        [betas ] = glmfit(X,y,'binomial','link','probit');
        betamat_misbnd(isub,it,:) = betas;

        d_diff = d2-d;
        d_diff = d_diff-min(d_diff);
        d_diff = d_diff./max(d_diff);

        y = ac(itrl1)>0;
        X = [td congresp  d d2 bt bl];
        %X = [td  d d2 bt bl];
        %X = [d ];
        X = X(itrl1,:);
        [betas ] = glmfit(X,y,'binomial','link','probit');
        betamat(isub,it,:) = betas;

        y = rt(itrl1);
        X = [td congresp d d2 bt bl];
        %X = [td  d d2 bt bl];
        %X = [d ];
        X = X(itrl1,:);
        [betas ] = glmfit(X,y,'normal','link','identity');
        betamat_rt(isub,it,:) = betas;
        
    end
    end
    
    
end
fprintf(' - done!\n');
%% averages measures for participants in both study 1 and 2
misbnd_acs = Switch_MergeExperiments(misbnd_ac,subpairs,explist,false,+Inf);
misbnd_resps = Switch_MergeExperiments(misbnd_resp,subpairs,explist,false,+Inf);
misbnd_resps_hat = Switch_MergeExperiments(misbnd_resp_hat,subpairs,explist,false,+Inf);
misbnd_rts = Switch_MergeExperiments(misbnd_rt,subpairs,explist,false,+Inf);

bmat_misbnd = Switch_MergeExperiments(betamat_misbnd,subpairs,explist,false,+Inf);
bmat_rt = Switch_MergeExperiments(betamat_rt,subpairs,explist,false,+Inf);
bmat = Switch_MergeExperiments(betamat,subpairs,explist,false,+Inf);
rout = Switch_MergeExperiments(outdata,subpairs,explist,false,+Inf);
mles = Switch_MergeExperiments(mle_vals,subpairs,explist,false,+Inf);

mles_misbind = Switch_MergeExperiments(mle_vals_misbind,subpairs,explist,false,+Inf);
tacs = Switch_MergeExperiments(tac,subpairs,explist,false,+Inf);
trts = Switch_MergeExperiments(trt,subpairs,explist,false,+Inf);
repfx = Switch_MergeExperiments(repeffect,subpairs,explist,false,+Inf);
outdat = Switch_MergeExperiments(outdata,subpairs,explist,false,+Inf);
dvdat = Switch_MergeExperiments(dvdata,subpairs,explist,false,+Inf);

mles  = squeeze(mles(:,:,:,1:2));
mles_misbind  = squeeze(mles_misbind(:,:,:,1:3));
%% plot accuracy by repeat
close all
clc
SetupPlotting();
figure('Position',[50 50 700 900]);

%-----------------%       
% plot acc and RT %
%-----------------%     
plotcol = [0.2000 0.3500 0.8000;  %blue (switch)
           0.3000 0.3000 0.3000;  %grey (repeat 1)
           0.7000 0.7000 0.7000]; %grey (repeat 2+)
measlabels = {'% correct' 'RT (s)'};       
for imeas = 1:2      
    subplot(5,2,[1 3]+imeas-1)
    if imeas == 1
        ylims = [0.5 1.0];
    else
        ylims = [0 1.5]; end
    xlims = [0.5 3.5];
    hold on
    plotdat = repfx(:,:,imeas);

    nt = size(plotdat,2);
    nx = size(plotdat,1);

    mubar = mean(plotdat,1);
    sdbar = std(bsxfun(@minus,plotdat,mean(plotdat,2)),[],1)./sqrt(size(plotdat,1));

    sdbar = std(plotdat,[],1)./sqrt(size(plotdat,1));
    for ibar = 1:nt
        boxchart(ibar*ones(nx,1),plotdat(:,ibar),'BoxFaceColor',plotcol(ibar,:),'Notch','off','MarkerColor',[0 0 0]);
        plot(ibar+[-0.25 +0.25],ones(1,2)*mubar(1,ibar),'k-','linewidth',0.5);
        for isub = 1:size(plotdat,1)
           hp = scatter(ibar+(rand(1,1)-0.5)*0.2,plotdat(isub,ibar),25,'markeredgecolor',[0 0 0],'markerfacecolor',plotcol(ibar,:),'MarkerFaceAlpha',0.5);
        end
    end
    offsetAxes()
    maxdat = max(plotdat(:))+range(plotdat(:))*0.05;
    ylabel(measlabels{imeas})
    xlabel('Trial Type')
    set(gca,'xtick',1:nt,'xticklabel',{'Switch' 'Repeat 1' 'Repeat 2+'});
    if imeas == 1
        set(gca, 'ytick',[0.5:.1:1],'yticklabel',{ '50'  '60'  '70'  '80'  '90' '100'}); end
    ylim(ylims)
    xlim(xlims)
    
    mubar = mubar*10^(1+imeas); 
    sdbar = sdbar*10^(1+imeas);
    fprintf('\n%s:\n',measlabels{imeas})
    fprintf('\nSwitch:%3.1f+%3.1f(s.e.m.)',mubar(1),sdbar(1))
    fprintf('\nRep1:%3.1f+%3.1f(s.e.m.)',mubar(2),sdbar(2))
    fprintf('\nRep2+:%3.1f+%3.1f(s.e.m.)',mubar(3),sdbar(3))
    fprintf('\n')
    plot_sig = 1;
    if plot_sig, 
        [~, pval,~,stat] = ttest(diff(plotdat(:,[1 2]),[],2));
        d = computeCohen_d(plotdat(:,1),plotdat(:,2),'paired');
        mysigstar(gca, [1.1 1.9], maxdat, pval,'k',0); 
        fprintf('\nSwitch vs. Rep1: t(%d)=%1.4f, p=%g, Cohen''s d=%1.2f',stat.df,stat.tstat,pval,d);
         [~, pval,~,stat] = ttest(diff(plotdat(:,[3 2]),[],2));
         d = computeCohen_d(plotdat(:,2),plotdat(:,3),'paired');
        fprintf('\nRep1 vs Rep2+: t(%d)=%1.4f, p=%g, Cohen''s d=%1.2f',stat.df,stat.tstat,pval,d);
        mysigstar(gca, [2.1 2.9], maxdat, pval,'k',0); 
        
        [~, pval,~,stat] = ttest(diff(plotdat(:,[1 3]),[],2));
        mysigstar(gca, [1.1 2.9], maxdat+range(plotdat(:))*0.20, pval,'k',0); 
         d = computeCohen_d(plotdat(:,1),plotdat(:,3),'paired');
        fprintf('\nSwitch vs Rep2+: t(%d)=%1.4f, p=%g, Cohen''s d=%1.2f',stat.df,stat.tstat,pval,d);
    end
    fprintf('\n\n')
    
end


%------------------------%
%plot psychometric curve %
%------------------------%

subplot(5,2,[7 9 ])

plotcol = [0.2000 0.3500 0.8000;  %blue (switch)
           0.5000 0.5000 0.5000]; %grey (repeat)

for icond = 1:2
    hold on
    %plot model fit
    m = squeeze(mles(:,icond,1:2));
    hold on
    plot([0 45],1-ones(1,2)*mean(m(:,2))/2,'k--','color',plotcol(icond,:),'linewidth',1)
    x = linspace(0,45,100);
    y = [];
    for isub = 1:size(m,1)
       y(isub,:) = (1-m(isub,2)).*normcdf(x,0,m(isub,1)) + m(isub,2)/2; 
    end
    cfg = [];
    cfg.alpha = 0.25;
    cfg.linewidth = 0.5;
    hp(icond) = plotpatch(y,x,plotcol(icond,:),cfg);
    %plot data
    cfg = [];
    cfg.dotsize = 25;
    plotDotsWithErrors(dvdat(:,:,icond),udv',plotcol(icond,:),cfg);
end
legend(hp,{'Switch' 'Repeat'},'location','SouthEast')
legend boxoff
set(gca,'xtick',round(udv))
set(gca, 'ytick',[0.0:.1:1],'yticklabel',{'0'  '10'  '20'  '30'  ...
        '40'  '50'  '60'  '70'  '80'  '90' '100'});
ylim([0.5 1])
ylabel('% correct')
xlabel('Probe Offset (\circ)')

%descriptives and stats for model fit %
mumle = mean(mles,1);
sdmle = std(mles,[],1)./sqrt(size(mles,1));
fprintf('\nModel Fit:\n')
fprintf('Switch: sigma=%2.1f+%2.1f, guess rate=%2.1f+%2.1f\n',mumle(1,1,1),sdmle(1,1,1),mumle(1,1,2)*100/2,sdmle(1,1,2)*100/2);
fprintf('Repeat: sigma=%2.1f+%2.1f, guess rate=%2.1f+%2.1f\n',mumle(1,2,1),sdmle(1,2,1),mumle(1,2,2)*100/2,sdmle(1,2,2)*100/2);

[h p c s] = ttest(squeeze(-diff(mles(:,:,1),[],2)));
d = computeCohen_d(mles(:,1,1),mles(:,2,1),'paired');
fprintf('\nSwitch vs Repeat sigma: t(%d)=%1.2f, p=%g, d=%1.2f',s.df,s.tstat,p,d);
[h p c s] = ttest(squeeze(-diff(mles(:,:,2),[],2)));
d = computeCohen_d(mles(:,1,2),mles(:,2,2),'paired');
fprintf('\nSwitch vs Repeat guess rate: t(%d)=%1.2f, p=%g, d=%1.2f',s.df,s.tstat,p,d);


%----------------------------------------------%
% plot effect of template distance on accuracy %
%----------------------------------------------%
subplot(5,2,8)
hold on
hpatch = [];
for icond = 1:2
    betas = permute(massGLM(utd,permute(mean(outdat(:,:,icond,1,:),5),[2 1 3 4])),[2 1]);
    cfg = [];
    cfg.linewidth = 0.5;
    cfg.alpha = 0.25;
    hpatch(icond)=plotpatch(betas(:,1,:,:) + betas(:,2,:,:)*utd',utd',plotcol(icond,:),cfg);
    cfg = [];
    cfg.dotsize = 25;
    plotDotsWithErrors(outdat(:,:,icond,1),utd',plotcol(icond,:),cfg);
end
legend(hpatch,{'Switch' 'Repeat'},'location','SouthWest')
legend boxoff
offsetAxes()
ylims = [0.70 0.85];
set(gca, 'ytick',[.70:.05:.85],'yticklabel',{'70' '75'  '80'  '85' });
%set(gca, 'ytick',[.35:.05:.95],'yticklabel',{'35' '' '45' '' '55' '' '65' '70' '75'  '80'  '85' '' '95'});
set(gca,'xtick',round(utd))
ylabel('% correct')
ylim(ylims)
xlim([0 90])

betas = permute(massGLM(utd,permute(outdat(:,:,:,1),[2 1 3])),[2 1 3 ]);
mubeta = mean(squeeze(betas(:,2,:)),1);
sdbeta = std(squeeze(betas(:,2,:)),[],1)./sqrt(size(betas,1));
[h p c s] = ttest(squeeze(betas(:,2,:)));
d1 = computeCohen_d(betas(:,2,1),betas(:,2,1)*0,'paired');
d2 = computeCohen_d(betas(:,2,2),betas(:,2,2)*0,'paired');

fprintf('\n\nEffect of template distance (Accuracy):')
fprintf('\nSwitch: %1.2f+%1.2f%% per 10 deg, t(%d)=%1.2f, p=%g, d=%1.2f',mubeta(1)*1000,sdbeta(1)*1000,s.df(1),s.tstat(1),p(1),d1)
fprintf('\nRepeat: %1.2f+%1.2f%% per 10 deg, t(%d)=%1.2f, p=%g, d=%1.2f',mubeta(2)*1000,sdbeta(2)*1000,s.df(2),s.tstat(2),p(2),d2)

[h p c s] = ttest(-diff(betas(:,2,:),[],3));
d = computeCohen_d(betas(:,2,1),betas(:,2,2),'paired');

fprintf('\nSwitch vs Repeat: t(%d)=%1.2f, p=%g,d=%1.2f',s.df(1),s.tstat(1),p(1),d);

[h p c s] = ttest(-diff(betas(:,1,:),[],3));
d = computeCohen_d(betas(:,1,1),betas(:,1,2),'paired');

fprintf('\nSwitch vs Repeat mean: t(%d)=%1.2f, p=%g,d=%1.2f',s.df(1),s.tstat(1),p(1),d);

%----------------------------------------%
% plot effect of template distance on RT %
%----------------------------------------%
subplot(5,2,10)
hold on
for icond = 1:2
    betas = permute(massGLM(utd,permute(outdat(:,:,icond,2),[2 1 3 4])),[2 1]);
    cfg = [];
    cfg.linewidth = 0.5;
    cfg.alpha = 0.25;
    plotpatch(betas(:,1,:,:) + betas(:,2,:,:)*utd',utd',plotcol(icond,:),cfg);
    cfg = [];
    cfg.dotsize = 25;
    plotDotsWithErrors(outdat(:,:,icond,2),utd',plotcol(icond,:),cfg);
end
offsetAxes()
ylims = [0.55 0.85];
ylim(ylims)
set(gca,'ytick',[0.6:.1:0.8]);
xlim([0 90])
set(gca,'xtick',round(utd))
ylabel('RT (s)')
xlabel('Template Distance (\circ)')


betas = permute(massGLM(utd,permute(outdat(:,:,:,2),[2 1 3])),[2 1 3 ]);
mubeta = mean(squeeze(betas(:,2,:)),1);
sdbeta = std(squeeze(betas(:,2,:)),[],1)./sqrt(size(betas,1));
[h p c s] = ttest(squeeze(betas(:,2,:)));
d1 = computeCohen_d(betas(:,2,1),betas(:,2,1)*0,'paired');
d2 = computeCohen_d(betas(:,2,2),betas(:,2,2)*0,'paired');

fprintf('\n\nEffect of template distance (RT):')
fprintf('\nSwitch: %1.2f+%1.2fms per 10 deg, t(%d)=%1.2f, p=%g, d=%1.2f',mubeta(1)*10000,sdbeta(1)*10000,s.df(1),s.tstat(1),p(1),d1)
fprintf('\nRepeat: %1.2f+%1.2fms per 10 deg, t(%d)=%1.2f, p=%g, d=%1.2f',mubeta(2)*10000,sdbeta(2)*10000,s.df(2),s.tstat(2),p(2),d2)

[h p c s] = ttest(-diff(betas(:,2,:),[],3));
d = computeCohen_d(betas(:,2,1),betas(:,2,2),'paired');

fprintf('\nSwitch vs Repeat slope: t(%d)=%1.2f, p=%g,d=%1.2f',s.df(1),s.tstat(1),p(1),d);

[h p c s] = ttest(-diff(betas(:,1,:),[],3));
d = computeCohen_d(betas(:,1,1),betas(:,1,2),'paired');

fprintf('\nSwitch vs Repeat mean: t(%d)=%1.2f, p=%g,d=%1.2f',s.df(1),s.tstat(1),p(1),d);
fprintf('\n')

%compare for each template distance:
[h p c s] = ttest(squeeze(diff(outdat,[],3)));
fprintf('\n\nAccuracy:')
for itd = 1:ntd
    d = computeCohen_d(outdat(:,itd,1,1),outdat(:,itd,2,1),'paired');
    fprintf('\n%d deg:  t(%d)=%1.2f, p=%g,d=%1.2f',round(utd(itd)),s.df(1,itd,1),s.tstat(1,itd,1),p(1,itd,1),d);
end
fprintf('\n\nRT:')
for itd = 1:ntd
    d = computeCohen_d(outdat(:,itd,1,2),outdat(:,itd,2,2),'paired');
    fprintf('\n%d deg:  t(%d)=%1.2f, p=%g,d=%1.2f',round(utd(itd)),s.df(1,itd,2),s.tstat(1,itd,2),p(1,itd,2),d);
end
fprintf('\n')

do_print = 0;
if do_print
   figpath =  [studypath '/figures/Figure1/']; if ~exist(figpath,'dir'), mkdir(figpath); end
   fname = sprintf('%s/Accuracy-RT-Psychometric-TemplateDistance',figpath);
   fname = sprintf('%s_%s',fname,datestr(now,'yyyymmdd-HHMM'));
   print(gcf,'-dpng',fname);
   print(gcf,'-depsc2',fname,'-painters');
   print(gcf,'-dpdf',fname,'-painters');
end
%% HELPER FUNCTIONS

function cost = costaccfun(P, design, data, verbose)
    % COSTACCFUN - Cost function for 2-parameter psychometric fit
    % P = [sigma, lapse_rate]
    sigma = P(1);
    lapse = P(2);
    
    % Predicted probability of clockwise response
    pred = (1 - lapse) .* normcdf(design, 0, sigma) + lapse/2;
    
    % Negative log likelihood
    cost = -sum(log(pred(data == 1))) - sum(log(1 - pred(data == 0)));
    
    if verbose
        fprintf('sigma=%.2f, lapse=%.3f, cost=%.2f\n', sigma, lapse, cost);
    end
end

function cost = costaccfun_misbind(P, design, data, verbose)
    % COSTACCFUN_MISBIND - Cost function for 3-parameter psychometric fit with swap errors
    % P = [sigma, lapse_rate, swap_rate]
    sigma = P(1);
    lapse = P(2);
    swap = P(3);
    
    % Predicted probability accounting for both items
    pred = (1 - lapse/2 - swap) .* normcdf(design(:, 1), 0, sigma) + ...
           swap .* normcdf(design(:, 2), 0, sigma) + lapse/2;
    
    % Negative log likelihood
    cost = -sum(log(pred(data == 1))) - sum(log(1 - pred(data == 0)));
    
    if verbose
        fprintf('sigma=%.2f, lapse=%.3f, swap=%.3f, cost=%.2f\n', ...
                sigma, lapse, swap, cost);
    end
end