classdef standTrialAnalysis < handle 
    % Preprocessing for standing balance data (dos Santos et. al. 2017)
    % Filename:	standTrialAnalysis.m
    % Authors:  Rika Sugimoto Dimitrova, James Hermus, Federico Tessari
    % Date:     27 Dec 2022
    % Description: When initalized, this class imports quiet-standing
    % subject data and preprocesses them.
    
    properties
        
        % Subject Specific Parameters
        subject % Subject number in study
        trial % Trial Number
        condition % Experimental condition {'OR','CR','OF','CF'}:
                  % OR: open eyes, rigid ground
                  % CR: closed eyes, rigid ground
                  % OF: open eyes, standing on foam
                  % CF: closed eyes, standing on foam
        condNames % Label for experimental conditions {'OR','CR','OF','CF'}
        
        sfrq % Sampling frequency (Hz)
        trialdur % Trial duration (s)
        
        % Measures
        t % Time (s) [Nx1]
        X % COM position (m): [Nx3] in cartesian coordinates (x,y,z)

    end
    
    methods
        function  [this] = standTrialAnalysis(subject, condition, trial)
            this.subject = subject;
            this.condition = condition;
            this.trial = trial;
            
            % During construction always run
            this.sfrq = 100; % Enforced sampling frequency (Hz)
            this.trialdur = 60; % Trial duration (s)
            this.condNames = {'OR','CR','OF','CF'};
            this.get_subjectData();
            
        end
        
        function [] = get_subjectData(this)
            warning('off','MATLAB:table:ModifiedAndSavedVarnames');
            % Load data:            
            tbl = (readtable(['BrownData/',...
                'dosSantos_data_2017/mkr/PDS' num2str(this.subject,'%02i')...
                this.condNames{this.condition} num2str(this.trial,'%i')...
                'mkr.txt']));
            warning('on','MATLAB:table:ModifiedAndSavedVarnames');
            com_PSD = table2array(tbl(:,{'COG_X','COG_Y','COG_Z'}));
            
            % dos Santos' coordinate frame:
            % x: AP, Front (+) - Hind (-)
            % y: vertical, Up (+) - Down (-)
            % z: ML, Right (+) - Left (-)
            
            com_pos(:,1) = (-1)*com_PSD(:,3);
            com_pos(:,2) = (-1)*com_PSD(:,1);
            com_pos(:,3) = com_PSD(:,2);
            
            dt = 1/this.sfrq;
            this.t = (dt:dt:this.trialdur)';
            this.X = com_pos-com_pos(1,:);
            
        end

    end
end