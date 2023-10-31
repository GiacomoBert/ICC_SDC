%% Compute ICC and SDC
% Example script to compute ICC and SDC at a channel level, for each time
% point, for EEG signal collected in two different condtions (like two time
% points for the same group of subjects). Should be adapted before re-use.

%% Author's info
% Author: Giacomo Bertazzoli
% Affiliations: IRCCS Istituto Centro San Giovanni di Dio Fatebenefratelli,Brescia IT; Centre for Mind/Brain Sciences, Rovereto IT;
% Date: May 2022
% e-mail: giacomo.bertazzoli@cognitiveneuroscience.it

%% Packages
% Requires:
% ICC https://www.mathworks.com/matlabcentral/fileexchange/21501-intraclass-correlation-coefficients

%% Data
% example data can be found here:
% https://gin.g-node.org/CIMeC/TMS-EEG_brain_connectivity_BIDS/src/master 

%%
clear
close all
clc
tic

%% load packages
addpath('/ICC'); % from https://www.mathworks.com/matlabcentral/fileexchange/21501-intraclass-correlation-coefficients

%% data parietal SOUND T1
addpath('/data_folder_T1'); %path to data in the 2 conditions
addpath('/data_folder_T2');

load('data_T1.mat'); %load data in the 2 conditions. In my case data were in fieldtrip format (keepsubj). Basically a 3D matrix, subjects X channels X time
load('data_T2.mat');

%% savepath
savepath=cd;

%% calculate ICC for each time point and electrode
% [this part takes long, but once computed and saved you can just comment it out]

%initialize output variables
ICC_calc11=zeros(63,600); %ICC calculated by me
ICC_calc13=zeros(63,600);

ICC1s=zeros(63,600); %ICC calculated using the ICC function 
ICC2s=zeros(63,600);
ICC3s=zeros(63,600);
ICC1k=zeros(63,600);
ICC2k=zeros(63,600);
ICC3k=zeros(63,600);

%
SDC_calc1=zeros(63,600); % SDC calculation
SDC_calc2=zeros(63,600);


for elec_count=1:size(data_T1.individual,2) %per electrode
    parfor time_count=1:size(data_T1.individual,3)%per instant in time
        
        %list instant X from all subjects in electrod Y in T0
        
        T0=data_T1.individual(:,elec_count,time_count);
        T1=data_T2.individual(:,elec_count,time_count);
        
        % compute anova
        s=length(T0);%numb of subjects
        N=numel([T0 T1]);%total numb of entries
        a=2; %number of levels
        
        %calculate DFs
        df_betw= a-1;
        df_with=N-a;
        df_subj=s-1;
        df_error=df_with-df_subj;
        df_total=N-1;
        
        % calculate SS
        SS_betw_trials= ((sum(T0)^2 + sum(T1)^2) /s) - (sum([T0; T1])^2/N);
        SS_with=sum([T0; T1].^2)-((sum(T0)^2 + sum(T1)^2)/s);
        SS_subj_betw= (sum((T0 + T1).^2)/a)- (sum([T0; T1])^2/N);
        SS_error=SS_with-SS_subj_betw;
        SS_total=SS_betw_trials+SS_with;
        
        %calculate MS
        MS_betw_trials=SS_betw_trials/df_betw;
        MS_with=SS_with/df_with;
        MS_subj_betw=SS_subj_betw/df_subj;
        MS_error=SS_error/df_error;
        MS_total=SS_total/df_total;
        
        % F (unnecessary)
        F=MS_betw_trials/MS_error;
        
        % ICC(1,1) and ICC(1,3) (from myself)
        SS_within_subjects=SS_error+SS_betw_trials;
        MS_within_subjects=SS_within_subjects/(df_error+df_betw);
        
        ICC_calc11(elec_count,time_count)=(MS_subj_betw-MS_within_subjects)/(MS_subj_betw+(a-1)*MS_within_subjects);
        ICC_calc13(elec_count,time_count)=(MS_subj_betw-MS_within_subjects)/(MS_subj_betw);
        
        %ICC test
        ICC1s(elec_count,time_count)=ICC(1,'single', [T0 T1]);
        ICC2s(elec_count,time_count)=ICC(2,'single', [T0 T1]);
        ICC3s(elec_count,time_count)=ICC(3,'single', [T0 T1]);
        ICC1k(elec_count,time_count)=ICC(1,'k', [T0 T1]);
        ICC2k(elec_count,time_count)=ICC(2,'k', [T0 T1]);
        ICC3k(elec_count,time_count)=ICC(3,'k', [T0 T1]);
        
        %SDC
        SDC_calc1(elec_count,time_count)=sqrt(MS_error) * 1.96 * sqrt(2);
        SDC_calc2(elec_count,time_count)=sqrt(MS_within_subjects) * 1.96 * sqrt(2);
        
    end
end
toc

%save output
save([savepath '\ICCs_SDCs'], 'ICC_calc11', 'ICC_calc13', 'ICC1s', 'ICC2s', ...
    'ICC3s', 'ICC1k', 'ICC2k', 'ICC3k', 'SDC_calc1', 'SDC_calc2');


%% plot ICC and channels
load([cd '\ICCs_SDCs.mat']); %load the results of the previous part 
load('data_T1_avg.mat'); %load the same data, but averaged across subjects channels X time
load('data_T2_avg.mat');

parfor elec_count=1:size(data_T1_avg.avg,1) %per electrode
    x=data_T1_avg.time; %time axes
    y1=[data_T1_avg.avg(elec_count,:); data_T2_avg.avg(elec_count,:) ]; %data axes 
    y2=ICC_calc13(elec_count,:); %ICC axes
    
    figure('Renderer', 'painters', 'Position', [10 10 1200 600]); %open a bigger window
    
    
    yyaxis left
    plot(x,y1,'k', 'LineWidth',2)
    set(gca,'YColor','k')
    
    hold on
    
    ylabel('Amplitude (\muV)');
    % yticks([-12 12]);
    % yticklabels({'-12', '12'});
    
    
    timeRect = [-0.001 0.006];
    xlim([-0.05 0.35]);
    ylim([-14 15]);
    ymax = 80;
    ymin = -40;
    rectangle('Position',[timeRect(1) ymin (timeRect(2)-timeRect(1)) ymax],'FaceColor',[0.7 0.7 0.7 0.3], 'LineStyle','none');
    
    
    
    yyaxis right
    plot(x, y2,'r', 'LineWidth',2)
    set(gca,'YColor','r')
    
    ylabel('ICC');
    ylim([-1 1]);
    
    xlabel('Time (ms)');
    set(gca,'box','off') %remove ticks on the figure
    xticklabels({'-50', '0', '', '100', '', '200', '', '300', ''})
    
    ax = gca;
    c = ax.Color;
    ax.FontWeight  = 'bold';
    ax.FontSize    = 30;
    set(gca,'TickDir','out'); % The only other option is 'in'
    
    title(['Channel: ' data_T1_avg.label{elec_count} ])
    hold off
    
    % save
    savefig([savepath '\' data_T1_avg.label{elec_count} ]);
    saveas(gcf,...
        [savepath '\' data_T1_avg.label{elec_count} ],...
        'tif');
    close
    
    
    
end

%% plot SDC and channels
load([cd '\ICCs_SDCs.mat']);
load('data_T1_avg.mat');
load('data_T2_avg.mat');

parfor elec_count=1:size(data_T1_avg.avg,1) %per electrode
    x=data_T1_avg.time;
    y1=[data_T1_avg.avg(elec_count,:); data_T2_avg.avg(elec_count,:) ];
    y2=SDC_calc1(elec_count,:)./...
        mean([data_T1_avg.avg(elec_count,:); ...
        data_T2_avg.avg(elec_count,:) ],1); %normalized SDC
    y3=SDC_calc1(elec_count,:);
    
    figure('Renderer', 'painters', 'Position', [10 10 1200 600]); %open a bigger window
    
    
    yyaxis left
    plot(x,y1,'k', 'LineWidth',2)
    set(gca,'YColor','k')
    hold on
    
    ylabel('Amplitude (\muV)');
    % yticks([-12 12]);
    % yticklabels({'-12', '12'});
    
    
    timeRect = [-0.001 0.006];
    xlim([-0.05 0.35]);
    ylim([-14 15]);
    ymax = 80;
    ymin = -40;
    rectangle('Position',[timeRect(1) ymin (timeRect(2)-timeRect(1)) ymax],'FaceColor',[0.7 0.7 0.7 0.3], 'LineStyle','none');
    
    
    
    yyaxis right
    plot(x, y3, 'b','LineWidth',2)
    set(gca,'YColor','b')
    
    ylabel('SDC (\muV)');
    ylim([0 10]);
    
    xlabel('Time (ms)');
    set(gca,'box','off') %remove ticks on the figure
    xticklabels({'-50', '0', '', '100', '', '200', '', '300', ''})
    
    ax = gca;
    c = ax.Color;
    ax.FontWeight  = 'bold';
    ax.FontSize    = 30;
    set(gca,'TickDir','out'); % The only other option is 'in'
    
    title(['Channel: ' data_T1_avg.label{elec_count} ])
    hold off
    % save
    savefig([savepath '\' data_T1_avg.label{elec_count} ]);
    saveas(gcf,...
        [savepath '\' data_T1_avg.label{elec_count} ],...
        'tif');
    close
    
    
    
end