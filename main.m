classdef main < handle
    % Filename:	main
    % Authors:  Federico Tessari, James Hermus, Rika Sugimoto Dimitrova
    % Date:		18 November 2022
    % Last modified: 22 March 2023
    % Description:	main program for analysis
    
    methods(Static)

        % Set default figure properties
        function set_figure_properties()

            set(0, 'DefaultLineLineWidth', 1);
            set(groot,'defaultAxesFontSize',8);
            set(0,'defaultfigurecolor',[1 1 1]); % white figure background
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(groot,'defaultTextInterpreter','latex');       
                                    
        end

    end

    properties

        % Boolean variables to choose which data to analyze
        bool_analyzeCrank
        bool_analyzeHand
        bool_analyzeStand
        bool_analyzeSim
        bool_runSim

        % Data variables
        %-- crank turning --%
        crankSubjDex  % Indices for crank turning experiment subjects
        crankTrialDex % Indices for crank turning experiment trials
        crankSfrq % Sampling frequency for crank experiments (Hz)
        thcp_r % Ensemble of crank position [Subject x trial x N] (Slow speed, CCW direction)
        t_ct    % Time vector, for crank data (s)
        CI_ct %Confidence Interval Crank
        R2_PSDct % R squared for PSD fitting

        %-- hand posture --%
        handSubjDex  % Indices for hand posture experiment subjects
        handTrialDex % Indices for hand posture experiment trials
        handSfrq     % Sampling frequency for hand posture experiments (Hz)
        xe   % Ensemble of x-coord of hand positions (m) [Subject x trial x N]
        ye   % Ensemble of y-coord of hand positions (m) [Subject x trial x N]
        t_hp    % Time vector, for hand data with vision (s)
        CI_hp %Confidence Interval hand position (x,y)
        R2_PSDhp % R squared for PSD fitting

        %-- quiet standing --%
        standSubjDex  % Indices for quiet standing experiment subjects
        standTrialDex % Indices for quiet standing experiment trials
        standSfrq     % Sampling frequency for quiet standing experiments (Hz)
        xcom   % Ensemble of x-coord of COM positions (m) [Subject x trial x N]
        ycom   % Ensemble of y-coord of COM positions (m) [Subject x trial x N]
        zcom   % Ensemble of z-coord of COM positions (m) [Subject x trial x N]
        t_qs   % Time vector, for quiet standing data (s)  
        CI_qs %Confidence Interval CoM position (x,y)
        R2_PSDqs % R squared for PSD fitting

        %----model output---$
        simSubjDex   % Indices for simulation "subjects" (1 by default)
        simTrialDex  % Indices for simulation trials
        simSfrq      % Sampling frequency for simulations (Hz)
        pos_sim1     % Ensemble of simulated crank angular position (rad) [Subject x trial x N]
        pos_sim2     % Ensemble of simulated hand angular position (rad) [Subject x trial x N]
        pos_sim3     % Ensemble of simulated stance angular position (rad) [Subject x trial x N]  
        t_sim1       % Time vector for crank turning simulation (s)
        t_sim2       % Time vector for hand posture simulation (s)
        t_sim3       % Time vector for quiet standing simulation (s)

        % Indices for representative subject(s) for whom to plot data
        crankMainSubjDex
        handMainSubjDex
        standMainSubjDex
        simMainSubjDex

        % Plot colors
        crank_color
        hand_color 
        stand_color

        % Power spectral density (PSD) computation parameters
        welch_overlap_frac % Fraction of overlap in the data windows used for Welch's method

        % PSD stats (95% confidence intervals and R-squared)
        PSD_stats_crank
        PSD_stats_hand
        PSD_stats_stand

        % Variance R-squared
        R2_crank_var  % R-squared values for crank-turning experiments, for each subject
        R2_hand_var_X % R-squared values for hand-posture experiments over different time intervals [Subject x max_time]
        R2_hand_var_Y % R-squared values for hand-posture experiments over different time intervals [Subject x max_time]
        R2_tot_hand_var_X % R-squared values over total hand-posture trial duration, for each subject
        R2_tot_hand_var_Y % R-squared values over total hand-posture trial duration, for each subject
        
    end

    methods
        function this = main()

            clc, close all
            dbstop if error
            restoredefaultpath;
            addpath('crank/','hand/','stand/','model/');

            %% Check crank_data is downloaded
            if(7~=exist('BrownData', 'dir'))
                error('The BrownData folder is missing. Make sure to download BrownData folder from dropbox and add it to your local repository. Information is provided in the README.md.');
            end

            % SET PARAMS
            this.welch_overlap_frac = 0.5;
            this.crankSfrq = 200;

            this.bool_analyzeCrank = true;
            this.bool_analyzeHand = true;
            this.bool_analyzeStand = true;
            this.bool_analyzeSim = true;
            this.bool_runSim = false;
            
            % Select subject(s) of interest to plot data
            % e.g. to select subjects 2 and 4, enter [2,4]
            this.crankMainSubjDex = 10; % Max 10 for crank-turning data
            this.handMainSubjDex = 10;  % Max 10 for hand-posture data
            this.standMainSubjDex = 4;  % Max 26 for quiet-standing data
            this.simMainSubjDex = 1;    % Max 1 for simulation

            this.crank_color = [0 0.4470 0.7410];
            this.hand_color = [0.8500 0.3250 0.0980]; 
            this.stand_color = [0.4660 0.6740 0.1880];

            % GET PREPROCESSED DATA    
            
            %-- crank turning --%
            if this.bool_analyzeCrank
                disp('Getting crank-turning data');
                crankAnalysis = crankCrossTrialAnalysis();
                this.thcp_r = crankAnalysis.ensemble{1}.thcp_r;
                dt_c = 1/this.crankSfrq;
                this.t_ct = dt_c.*(1:size(this.thcp_r(1,1,:),3))';
                this.crankSubjDex = 1:size(this.thcp_r,1);
                this.crankTrialDex = 1:size(this.thcp_r,2);
            end%bool_analyzeCrank

            %-- hand posture --%
            if this.bool_analyzeHand
                disp('Getting hand-posture data');
                handAnalysis = handCrossTrialAnalysis();
                this.xe = handAnalysis.ensemble{1}.x;
                this.ye = handAnalysis.ensemble{1}.y; 
                this.t_hp = handAnalysis.time{1};  
                this.handSubjDex = 1:size(this.xe,1);
                this.handTrialDex = 1:size(this.xe,2);
                this.handSfrq = round(1/(this.t_hp(2)-this.t_hp(1)));
            end%bool_analyzeHand

            %-- quiet standing --%
            if this.bool_analyzeStand
                disp('Getting quiet-standing data');
                standAnalysis = standCrossTrialAnalysis();
                this.xcom = standAnalysis.ensemble{1}.xcom;
                this.ycom = standAnalysis.ensemble{1}.ycom; 
                this.zcom = standAnalysis.ensemble{1}.zcom; 
                this.t_qs  = standAnalysis.time{1};
                this.standSubjDex = 1:size(this.xcom,1);
                this.standTrialDex = 1:size(this.xcom,2);
                this.standSfrq = round(1/(this.t_qs(2)-this.t_qs(1)));
            end%bool_analyzeStand

            %-- simulations --%
            if this.bool_runSim
                N_trials = 200;% Number of iterations to run
                systems = 1:3; % Systems to simulate - 1: crank turning, 2: hand posture, 3: quiet standing
                datafiles = {'out_crank_new','out_hand_new','out_stand_new'}; % Simulation data files for the 3 systems
                %-- crank turning simulation --%
                if ismember(1,systems)
                    body_model = 1;
                    control_mode = 2;
                    disp('Running crank-turning simulations');
                    modelSim = simulateModel(body_model, control_mode, N_trials);
                    out = modelSim.simOut;
                    save(['./BrownData/modelling/' datafiles{1} '.mat'],'out');
                end
                %-- hand posture simulation --%
                if ismember(2,systems)
                    body_model = 1;
                    control_mode = 1;
                    disp('Running hand-posture simulations');
                    modelSim = simulateModel(body_model, control_mode, N_trials);
                    out = modelSim.simOut;
                    save(['./BrownData/modelling/' datafiles{2} '.mat'],'out');
                end
                %-- quiet standing simulation --%
                if ismember(3,systems)
                    body_model = 2;
                    control_mode = 1;
                    disp('Running quiet-standing simulations');
                    modelSim = simulateModel(body_model, control_mode, N_trials);
                    out = modelSim.simOut;
                    save(['./BrownData/modelling/' datafiles{3} '.mat'],'out');
                end
            end%bool_runSim

            if this.bool_analyzeSim
                if ~this.bool_runSim
                    datafiles = {'out_v_crank','out_v_hand','out_v_stand'};
                    N_trials = 200;
                    systems = 1:3;
                end
                disp('Getting simulation output');
                modelAnalysis = ...
                    modelCrossTrialAnalysis(systems,1:N_trials,datafiles);
                this.pos_sim1 = modelAnalysis.ensemble{1}.pos;
                this.pos_sim2 = modelAnalysis.ensemble{2}.pos; 
                this.pos_sim3 = modelAnalysis.ensemble{3}.pos; 
                this.t_sim1  = modelAnalysis.time{1};
                this.t_sim2  = modelAnalysis.time{2};
                this.t_sim3  = modelAnalysis.time{3};
                this.simSubjDex = 1:size(this.pos_sim1,1);
                this.simTrialDex = 1:size(this.pos_sim1,2);
                this.simSfrq = round(1/(this.t_sim1(2)-this.t_sim1(1)));
            end%bool_analyzeSim
            

            % PLOT PROCESSED DATA
            this.set_figure_properties();
            this.plot_main();
%             this.compute_PSD_stats();
            this.plot_supp();
%             this.compute_var_R2();

        end

        % Plots for paper, Main Text
        function plot_main(this)
            lineSpecs.width = 1;

            disp('Plotting main figures');
        
            %-- crank turning --%
            if this.bool_analyzeCrank                
            for subj = this.crankMainSubjDex % select one subject
                hf = figure('Units','centimeters','Position',[10 5 5.7 9]);
                subplot(311)
                plot(this.t_ct,squeeze(this.thcp_r(subj,:,:)))
                ylabel('$\theta$ [rad]')
                xlabel('time [s]')
                set(gca,'Box','off')
                grid on
                subplot(312)
                plot(this.t_ct,var(squeeze(this.thcp_r(subj,:,:))),'Color',this.crank_color)
                ylabel('variance [rad$^2$]')
                xlabel('time [s]')
                set(gca,'Box','off')
                grid on
                subplot(313)
                lineSpecs.color = this.crank_color;
                this.plot_ensemble_PSD(...
                    this.thcp_r, this.crankSfrq, this.t_ct(end),...
                    hf.Number, lineSpecs,subj, 'crank');
                ylim([-60,0])
                xlim([0.07 2])
                yticks(-60:20:0)
                ylabel('$PSD$ [dB/Hz]')
                xlabel('frequency [Hz]')
                sgtitle('Crank Turning')
            end
            end
            
            %-- hand posture --%
            if this.bool_analyzeHand
            for subj = this.handMainSubjDex % select one subject    
                hf = figure('Units','centimeters','Position',[10 5 5.7 9]);
                subplot(311)
                hold on
                plot(this.t_hp,squeeze(this.xe(subj,this.handTrialDex,:)))
                %plot(this.t_hp,squeeze(this.ye(subj,this.handTrialDex,:)))
                ylabel('hand trajectory [m]')
                set(gca,'Box','off')
                grid on
                subplot(312)
                hold on
                plot(this.t_hp,var(squeeze(this.xe(subj,:,:))),'Color',this.hand_color)
                %plot(this.t_hp,var(squeeze(this.ye(subj,:,:))),'Color',this.hand_color)
                ylabel('variance [m$^2$]')
                xlabel('time [s]')
                set(gca,'Box','off')
                grid on
                subplot(313)
                hold on
                lineSpecs.color = this.hand_color;
                this.plot_ensemble_PSD(...
                    this.xe, this.handSfrq, this.t_hp(end),...
                    hf.Number, lineSpecs,subj, 'hand');
                %this.plot_ensemble_PSD(...
                %   this.ye, this.handSfrq, this.t_hp(end),...
                %   hf.Number, lineSpecs, subj, 'hand');
                grid on
                xlim([0.01 1])
                ylim([-100,-40])
                yticks(-100:20:-40)
                ylabel('$PSD$ [dB/Hz]')
                xlabel('frequency [Hz]')
                sgtitle('Hand Posture')
            end
            end

            %-- quiet standing --%
            if this.bool_analyzeStand
            for subj = this.standMainSubjDex % select one subject
                hf = figure('Units','centimeters','Position',[10 5 5.7 6]);
                subplot(211)
                hold on
                plot(this.t_qs,squeeze(this.xcom(subj,this.standTrialDex,:)))
                %plot(this.t_qs,squeeze(this.ycom(subj,this.standTrialDex,:)))
                ylabel('$COM$ [m]')
                xlabel('time [s]')
                set(gca,'Box','off')
                grid on
                subplot(212)
                hold on
                lineSpecs.color = this.stand_color;
                this.plot_ensemble_PSD(...
                    this.xcom, this.standSfrq, this.t_qs(end),...
                    hf.Number,lineSpecs, subj, 'stand');
                %this.plot_ensemble_PSD(...
                %    this.ycom, this.standSfrq, this.t_qs(end),...
                %    hf.Number,lineSpecs, subj, 'stand');
                grid on
                xlim([0.03 0.4])
                ylim([-80,-40])
                yticks(-80:20:-40)
                ylabel('$PSD$ [dB/Hz]')
                xlabel('frequency [Hz]')
                sgtitle('Quiet Standing')
            end
            end

            %-- model simulation --%
            if this.bool_analyzeSim
                subj = this.simMainSubjDex;

                hf = figure('Units','centimeters','Position',[10 5 12.1 9]);

                % plot a subset of the trials
                N_trials2plot = 15;
                iter = length(this.simTrialDex);
                sub_iter = round(iter/N_trials2plot);
                trials2plot = 1:sub_iter:iter;

                %-- crank turning simulation --%
                wref = 5*2*pi/60; % reference velocity [rad/s];
                % trajectory
                subplot(331)
                hold on
                % align crank data after removing transients
                index_start = 4*401;
                t_plot = this.t_sim1(index_start:end) - this.t_sim1(index_start);
                pos_sim1_aligned = this.pos_sim1(subj,:,index_start:end)...
                    - this.pos_sim1(subj,:,index_start);
                plot(t_plot,...
                    squeeze(pos_sim1_aligned(subj,trials2plot,:))...
                    - wref*t_plot'.*ones(size(wref*t_plot')))
                plot(t_plot,(wref*t_plot - wref*t_plot),'k--')
                set(gca,'Box','off')
                grid on
                xlabel('time [s]')
                ylabel('position [rad]')
                title('Crank Turning')
                % variance
                subplot(334)
                plot(t_plot,...
                    var(squeeze(pos_sim1_aligned(subj,:,:))),...
                    'Color',this.crank_color)
                ylabel('variance [rad$^2$]')
                xlabel('time [s]');
                set(gca,'Box','off')
                grid on
                % PSD
                subplot(337);
                lineSpecs.color = this.crank_color;
                this.plot_ensemble_PSD(...
                    pos_sim1_aligned(subj,:,:), this.simSfrq, t_plot(end), ...
                    hf.Number, lineSpecs, subj, 'crank_sim');                
                xlim([0.01 1])
                ylim([-40 0])
                yticks(-40:20:0)
                ylabel('$PSD$ [dB/Hz]')
                xlabel('frequency [Hz]')            

                %-- hand posture simulation --%
                th_thresh = 2*pi/180; % threshold position error [rad]
                % trajectory
                subplot(332)
                hold on
                plot(this.t_sim2, squeeze(this.pos_sim2(subj,trials2plot,:)))
                plot(this.t_sim2,ones(size(this.t_sim2))*th_thresh,'k--')
                plot(this.t_sim2,ones(size(this.t_sim2))*(-1)*th_thresh,'k--')
                set(gca,'Box','off')
                grid on
                xlabel('time [s]')
                title('Hand Posture')
                % variance
                subplot(335)
                plot(this.t_sim2,...
                    var(squeeze(this.pos_sim2(subj,:,:))),...
                    'Color',this.hand_color)
                xlabel('time [s]')
                set(gca,'Box','off')
                grid on
                % PSD
                subplot(338)
                lineSpecs.color = this.hand_color;
                this.plot_ensemble_PSD(...
                    this.pos_sim2, this.simSfrq, ...
                    this.t_sim2(end),...
                    hf.Number, lineSpecs, subj, 'hand_sim');                
                xlim([0.01 1])
                ylim([-60 -20])
                yticks(-60:20:-20)
                xlabel('frequency [Hz]')

                %-- quiet standing simulation --%
                % trajectory
                subplot(333)
                hold on;
                plot(this.t_sim3, squeeze(this.pos_sim3(subj,trials2plot,:)))
                plot(this.t_sim3,ones(size(this.t_sim3))*th_thresh,'k--')
                plot(this.t_sim3,ones(size(this.t_sim3))*(-1)*th_thresh,'k--')
                set(gca,'Box','off')
                grid on
                xlabel('time [s]')
                title('Quiet Standing')
                % variance
                subplot(336)
                plot(this.t_sim2,...
                    var(squeeze(this.pos_sim3(subj,:,:))),...
                    'Color',this.stand_color)
                xlabel('time [s]')
                set(gca,'Box','off')
                grid on
                % PSD
                subplot(339)
                lineSpecs.color = this.stand_color;
                this.plot_ensemble_PSD(...
                    this.pos_sim3, this.simSfrq, ...
                    this.t_sim3(end),...
                    hf.Number, lineSpecs, subj, 'stand_sim');                
                xlim([0.01 1])
                ylim([-60 -20])
                yticks(-60:20:-20)
                xlabel('frequency [Hz]')
            end%bool_analyze_models

        end%plot_main

        % Plots for paper, Supplementary Materials
        % Uses stats computed from compute_PSD_stats,
        % saved in a spreadsheet under 'BrownData/PSD_stats.xlsx'
        function plot_supp(this)
            color_blue = [0 0.4470 0.7410];
            color_red = [0.8500 0.3250 0.0980];

            disp('Plotting supplementary confidence interval figures');

            %Confidence Interval Plots
            CI_crank = readtable('BrownData/PSD_stats.xlsx','Sheet','Crank');
            CI_hand = readtable('BrownData/PSD_stats.xlsx','Sheet','Hand');
            CI_stand = readtable('BrownData/PSD_stats.xlsx','Sheet','Stand');

            figure('Units','centimeters','Position',[10 5 5.7 5.7])
            N = length(CI_crank.CI_mean);
            errorbar(1:N,CI_crank.CI_mean,...
                CI_crank.CI_mean-CI_crank.CI_min,...
                CI_crank.CI_mean-CI_crank.CI_max,...
                'LineWidth',2,'Color',color_blue,'LineStyle','none'),
            hold on            
            plot(1:N,-20*ones(N,1),'--','Color',color_blue)
            xlim([0 N+1])
            xticks(1:1:N)
            xlabel('Subjects')
            ylabel('$PSD$ Slope CI [dB/dec]')
            ylim([-40 0])
            set(gca,'Box','off')            
            title('Crank Turning')

            figure('Units','centimeters','Position',[10 5 5.7 9])
            subplot(211)
            x_indices = find([CI_hand.direction{:}]=='x');
            N = length(x_indices);
            errorbar(1:N,CI_hand.CI_mean(x_indices),...
                CI_hand.CI_mean(x_indices)-CI_hand.CI_min(x_indices),...
                CI_hand.CI_mean(x_indices)-CI_hand.CI_max(x_indices),...
                'LineWidth',2,'Color',color_blue,'LineStyle','none'),
            hold on            
            plot(1:N,-20*ones(N,1),'--','Color',color_blue)
            xlim([0 N+1])
            xticks(1:1:N)
            ylim([-40 0])
            ylabel('$PSD_x$ Slope CI [dB/dec]')
            set(gca,'Box','off')    
            title('Hand Posture')
            
            subplot(212)
            y_indices = find([CI_hand.direction{:}]=='y');
            N = length(y_indices);
            errorbar(1:N,CI_hand.CI_mean(y_indices),...
                CI_hand.CI_mean(y_indices)-CI_hand.CI_min(y_indices),...
                CI_hand.CI_mean(y_indices)-CI_hand.CI_max(y_indices),...
                'LineWidth',2,'Color',color_red,'LineStyle','none'),
            hold on            
            plot(1:N,-20*ones(N,1),'--','Color',color_red)
            xlim([0 N+1])
            xticks(1:1:N)
            ylim([-40 0])
            xlabel('Subjects')
            ylabel('$PSD_y$ Slope CI [dB/dec]')
            set(gca,'Box','off')    


            figure('Units','centimeters','Position',[10 5 12.1 9])
            subplot(211)
            x_indices = find([CI_stand.direction{:}]=='x');
            N = length(x_indices);
            errorbar(1:N,CI_stand.CI_mean(x_indices),...
                CI_stand.CI_mean(x_indices)-CI_stand.CI_min(x_indices),...
                CI_stand.CI_mean(x_indices)-CI_stand.CI_max(x_indices),...
                'LineWidth',2,'Color',color_blue,'LineStyle','none'),
            hold on
            plot(1:N,-20*ones(N,1),'--','Color',color_blue)
            xlim([0 N+1])
            xticks(1:1:N)
            ylim([-40 0])
            ylabel('$PSD_x$ Slope CI [dB/dec]')
            set(gca,'Box','off')    
            title('Quiet Standing')
            
            subplot(212)
            y_indices = find([CI_stand.direction{:}]=='y');
            N = length(y_indices);
            errorbar(1:N,CI_stand.CI_mean(y_indices),...
                CI_stand.CI_mean(y_indices)-CI_stand.CI_min(y_indices),...
                CI_stand.CI_mean(y_indices)-CI_stand.CI_max(y_indices),...
                'LineWidth',2,'Color',color_red,'LineStyle','none'),
            hold on            
            plot(1:N,-20*ones(N,1),'--','Color',color_red)
            xlim([0 N+1])
            xticks(1:1:N)
            ylim([-40 0])
            xlabel('Subjects')
            ylabel('$PSD_y$ Slope CI [dB/dec]')  
            set(gca,'Box','off')
        end%plot_supp     

        % Compute statistics for best-fit PSD slopes, and save to file
        % (to be used by this.plot_supp to plot confidence intervals)
        function compute_PSD_stats(this)
            disp('Computing PSD stats');

            filename = 'BrownData/PSD_stats.xlsx';
            %-- crank turning --%
            if this.bool_analyzeCrank
            for subj = this.crankSubjDex
                [~,~,f_est,~,CI_slope,est_sl,R2_PSDo,R2_PSDa] ...
                    = this.fit_PSD_slope(...
                    this.thcp_r, this.crankSfrq, this.t_ct(end),...
                    subj, 'crank');
                subject(subj,1) = subj;
                freq_range_Hz{subj,1} = [num2str(round(f_est(1),1,'significant')) '-' num2str(round(f_est(end),1,'significant'))];
                CI_min(subj,1) = CI_slope(1);
                CI_max(subj,1) = CI_slope(2);
                CI_mean(subj,1) = est_sl;
                R2_ord(subj,1) = R2_PSDo;
                R2_adj(subj,1) = R2_PSDa;
            end
            this.PSD_stats_crank = table(subject,freq_range_Hz,CI_min,CI_max,CI_mean,R2_ord,R2_adj);
            writetable(this.PSD_stats_crank,filename,'Sheet','Crank','Range','A1');
            clearvars subject freq_range_Hz CI_min CI_max CI_mean R2_ord R2_adj
            end     

            %-- hand posture --%
            if this.bool_analyzeHand
            for subj = this.handSubjDex
                directions = ['x','y'];
                N_dir = length(directions);
                for i = 1:N_dir
                    switch directions(i)
                        case 'x'
                            [~,~,f_est,~,CI_slope,est_sl,R2_PSDo,R2_PSDa]...
                                = this.fit_PSD_slope(...
                                this.xe, this.handSfrq, this.t_hp(end),...
                                subj, 'hand');
                        case 'y'
                            [~,~,f_est,~,CI_slope,est_sl,R2_PSDo,R2_PSDa]...
                                = this.fit_PSD_slope(...
                                this.ye, this.handSfrq, this.t_hp(end),...
                                subj, 'hand');
                    end
                    subject(N_dir*(subj-1)+i,1) = subj;
                    direction(N_dir*(subj-1)+i,1) = directions(i);
                    freq_range_Hz{N_dir*(subj-1)+i,1} = [num2str(round(f_est(1),1,'significant')) '-' num2str(round(f_est(end),1,'significant'))];
                    CI_min(N_dir*(subj-1)+i,1) = CI_slope(1);
                    CI_max(N_dir*(subj-1)+i,1) = CI_slope(2);
                    CI_mean(N_dir*(subj-1)+i,1) = est_sl;
                    R2_ord(N_dir*(subj-1)+i,1) = R2_PSDo;
                    R2_adj(N_dir*(subj-1)+i,1) = R2_PSDa;   
                end
            end
            this.PSD_stats_hand = table(subject,direction,freq_range_Hz,CI_min,CI_max,CI_mean,R2_ord,R2_adj);
            writetable(this.PSD_stats_hand,filename,'Sheet','Hand','Range','A1');
            clearvars subject direction freq_range_Hz CI_min CI_max CI_mean R2_ord R2_adj
            end

            %-- quiet standing --%
            if this.bool_analyzeStand
            for subj = this.standSubjDex
                directions = ['x','y'];
                N_dir = length(directions);
                for i = 1:N_dir
                    switch directions(i)
                        case 'x'
                            [~,~,f_est,~,CI_slope,est_sl,R2_PSDo,R2_PSDa]...
                                = this.fit_PSD_slope(...
                                this.xcom, this.standSfrq, this.t_qs(end),...
                                subj, 'stand');
                        case 'y'
                            [~,~,f_est,~,CI_slope,est_sl,R2_PSDo,R2_PSDa]...
                                = this.fit_PSD_slope(...
                                this.ycom, this.standSfrq, this.t_qs(end),...
                                subj, 'stand');
                    end
                    subject(N_dir*(subj-1)+i,1) = subj;
                    direction(N_dir*(subj-1)+i,1) = directions(i);
                    freq_range_Hz{N_dir*(subj-1)+i,1} = [num2str(round(f_est(1),1,'significant')) '-' num2str(round(f_est(end),1,'significant'))];
                    CI_min(N_dir*(subj-1)+i,1) = CI_slope(1);
                    CI_max(N_dir*(subj-1)+i,1) = CI_slope(2);
                    CI_mean(N_dir*(subj-1)+i,1) = est_sl;
                    R2_ord(N_dir*(subj-1)+i,1) = R2_PSDo;
                    R2_adj(N_dir*(subj-1)+i,1) = R2_PSDa;   
                end
            end
            this.PSD_stats_stand = table(subject,direction,freq_range_Hz,CI_min,CI_max,CI_mean,R2_ord,R2_adj);
            writetable(this.PSD_stats_stand,filename,'Sheet','Stand','Range','A1');
            clearvars subject direction freq_range_Hz CI_min CI_max CI_mean R2_ord R2_adj
            end
            
        end

        % Compute R-squared stats for the variance plots
        function compute_var_R2(this)
            disp('Computing variance R-squared');

            %-- crank --%
            if this.bool_analyzeCrank
            for subj = this.crankSubjDex % select one subject
                lm = fitlm(this.t_ct,var(squeeze(this.thcp_r(subj,:,:))));
                this.R2_crank_var(subj,1) = lm.Rsquared.Ordinary;
            end
            end
            %-- hand posture --%
            if this.bool_analyzeHand
            t_min = 1;
            t_max = 200;
            for subj = this.handSubjDex % select one subject
                varX = var(squeeze(this.xe(subj,:,:)));
                varY = var(squeeze(this.ye(subj,:,:)));

                directions = ['x','y'];
                N_dir = length(directions);
                for i = 1:N_dir                    
                    lm_tot = fitlm(this.t_hp,varX);
                    this.R2_tot_hand_var_X(subj,1) = lm_tot.Rsquared.Ordinary;
                    time_s = 1;
                    for idx = (t_min:t_max)*this.handSfrq 
                        lm = fitlm(this.t_hp(1:idx),varX(1:idx));
                        this.R2_hand_var_X(subj,time_s) = lm.Rsquared.Ordinary;
                        time_s = time_s+1;
                    end
                    lm_tot = fitlm(this.t_hp,varY);
                    this.R2_tot_hand_var_Y(subj,1) = lm_tot.Rsquared.Ordinary;
                    time_s = 1;
                    for idx = (t_min:t_max)*this.handSfrq
                        lm = fitlm(this.t_hp(1:idx),varY(1:idx));
                        this.R2_hand_var_Y(subj,time_s) = lm.Rsquared.Ordinary;
                        time_s = time_s+1;
                    end
                end
            end
            end
        end

        % Plot the power spectral density (PSD) for the given ensemble
        function [mean_plot,slope_plot,f_est,px_est,CI_slope,est_sl,R2_PSDo,R2_PSDa] = ...
                plot_ensemble_PSD(this,ensemble,sfrq,T_max_s,ifig,lineSpecs,subj,datalabel)

            if nargin < 5 || isempty(ifig)
                figure;hold on;
            else
                figure(ifig);hold on;
            end
            
            if nargin < 6 || isempty(lineSpecs)
                lineSpecs.color = [0 0.4470 0.7410];
                lineSpecs.width = 1;
            end

            [p_mean,f_Hz,f_est,px_est,CI_slope,est_sl,R2_PSDo,R2_PSDa] = ...
                fit_PSD_slope(this,ensemble,sfrq,T_max_s,subj,datalabel);

            p_mean_log10 = 10*log10(p_mean);
            mean_plot = ...
                plot(f_Hz,p_mean_log10,'-',...
                'Linewidth',lineSpecs.width,...
                'Color',lineSpecs.color);

            if ( strcmp(datalabel,'crank_sim') ...
                || strcmp(datalabel,'hand_sim')...
                || strcmp(datalabel,'stand_sim') )
                % plot reference -20 dB/dec slope
                offset_dB = 3;
                slope_20dB = ...
                    -20*(log10(f_est)-log10(f_est(1))) + px_est(1) - offset_dB;
                slope_plot = ...
                    plot(f_est,slope_20dB,'k--','Linewidth',lineSpecs.width);
            else
                % plot best-fit line
                slope_plot = ...
                    plot(f_est,px_est,'k--','Linewidth',lineSpecs.width);
            end

            xlabel('Frequency (Hz)');set(gca,'XScale','log');
            xticks([0.01,0.1,1,10,100]);

            set(gca,'Box','off');
            grid on;
                                   
        end

        % Get power spectral density (PSD), best-fit line,
        % and related stats (95% confidence interval (CI) and R-squared)
        function [p_mean,f_Hz,f_est,px_est,CI_slope,est_sl,R2_PSDo,R2_PSDa] = ...
                fit_PSD_slope(this,ensemble,sfrq,T_max_s,subj,datalabel)

            [p_all,f_Hz] = get_PSD(this,ensemble,sfrq,T_max_s,subj,datalabel);    
            
            p_mean = mean(p_all,2);

            f_in = 0.01;
            f_out = 1;
            switch datalabel
                case 'crank'
                    f_in = 0.07;
                    f_out = 1;
                case 'hand'
                    if subj == 2 || subj == 8
                        f_in = 0.03;
                        f_out = 0.15;
                    elseif subj == 5
                        f_in = 0.04;
                        f_out = 0.2;
                    else
                        f_out = 0.1;
                    end
                case 'stand'
                    f_in = 0.033;
                    f_out = 0.4;
                case 'stand_sim'
                    f_out = 0.1;
            end
            
            [std_err,f_est,px_est,est_sl,R2_PSDo,R2_PSDa] = ...
                PSD_slope_reg(this,f_Hz,p_mean,f_in,f_out);
            Nci = length(px_est);
            CI_slope = conf_int_slope(this,est_sl,95,std_err,Nci);
            %disp(est_sl);
            %disp(CI_slope);
            %disp(R2_PSDo);
        end

        % Compute the PSD for a given subject in an ensemble
        function [p_all,f_Hz] = get_PSD(this,ensemble,sfrq,T_max_s,subj,datalabel)          

            %disp('Computing PSD');

            N_subj = size(ensemble,1);
            N_trials = size(ensemble,2);
            window_n = round(T_max_s*sfrq);

            p_all = zeros(window_n + 1 - 4,N_subj*N_trials);

            count = 1;                          
            for trial = 1:N_trials
                meas_temp = squeeze(ensemble(subj,trial,:));
                if strcmp(datalabel,'crank') || strcmp(datalabel,'crank_sim')
                    meas_detrend = meas_temp-mean(squeeze(ensemble(subj,:,:)))';
                else
                    meas_detrend = detrend(meas_temp,0);  
                end
            
                window_n = round(T_max_s*sfrq);
            
                [p,f_Hz] =...
                    cpsd(meas_detrend,meas_detrend,...
                    rectwin(window_n),window_n*this.welch_overlap_frac,window_n*2,sfrq);
                p_all(:,count) = p(5:end);
                f_Hz = f_Hz(5:end);
            
                count = count+1;
            end
        end

        % Get best-fit line and R-squared for a single low-frequency PSD
        function [std_err,f_est,px_est_dB,est_sl,R2_PSDo,R2_PSDa] = ...
                PSD_slope_reg(this,f,px,f_in,f_out)

            idx = find(f>=f_in);
            idx_in = idx(1);
            clearvars idx

            idx = find(f<=f_out);
            idx_end = idx(end);
            clearvars idx

            % Fit linear model
            mdl = fitlm(log10(f(idx_in:idx_end)),10*log10(px(idx_in:idx_end)));
            
            % Get best-fit line
            est_sl = mdl.Coefficients{2,1}; % estimated slope
            f_est = f(idx_in:idx_end);
            px_est_dB = est_sl*log10(f(idx_in:idx_end))+mdl.Coefficients{1,1};
            R2_PSDo = mdl.Rsquared.Ordinary;
            R2_PSDa = mdl.Rsquared.Adjusted;

            px_dB = 10*log10(px(idx_in:idx_end));
            f_log10 = log10(f(idx_in:idx_end));
            f_mean_log10 = mean(f_log10);
            N = length(f(idx_in:idx_end));

            % Standard Error
            std_err = sqrt( (1/(N-2))*sum((px_dB-px_est_dB).^2)./sum((f_log10-f_mean_log10).^2) );
        end

        function CI_slope = conf_int_slope(this,slope,interval,st_er,N)

            range = (100-interval)/200;

            ts = tinv([range  1-range],N-2);      % T-Score
            CI_slope = slope + ts*st_er;          % Confidence Intervals
        end

    end

                            
end

