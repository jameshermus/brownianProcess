classdef handTrialAnalysis < handle 
    % Data preprocessing for hand posture control experiments
    % Filename:	handTrialAnalysis.m
    % Authors:  Rika Sugimoto Dimitrova, Federico Tessari, James Hermus
    % Date:     19 Dec 2022
    % Description: When initalized, this class imports subject data and
    % preprocesses them.
    
    properties
        
        % Subject Specific Parameters
        subject % Subject number in study
        trial % Trial Number
        condition % Experimental condition {1,2}: 1: with visual feedback (v), 2: no visual feedback (nv)
        condNames % Label for experimental conditions {'v', 'nv'}
        
        sfrq % Sampling frequency (Hz)
        trialdur % Trial duration (s)
        
        % Measures
        t % Time (s) [Nx1]
        X % Handle position (m): [Nx2] in cartesian crank coordinates

    end
    
    methods
        function  [this] = handTrialAnalysis(subject, condition, trial)
            this.subject = subject;
            this.condition = condition;
            this.trial = trial;
            
            % During construction always run
            this.sfrq = 200; % Enforced sampling frequency (Hz)
            this.trialdur = 4*60; % Trial duration (s)
            this.condNames = {'v','nv'};
            this.get_subjectData();
            
        end
        
        function [] = get_subjectData(this)        
            % Load data:
            data = load(['BrownData/hand_post/',...
                'Subject' int2str(this.subject) '/' ...
                'subj' int2str(this.subject) '_trial' int2str(this.trial),'.mat']);%'_' this.condNames{this.condition} ,'.mat']);
            
            % Resample
            xe_rs = resample(data.posX_m,data.time,this.sfrq);
            ye_rs = resample(data.posY_m,data.time,this.sfrq);
            trs = (0:1/this.sfrq:1/this.sfrq*length(xe_rs)-1/this.sfrq)';
            
            N = length(trs);
            end_t = trs(N);
            init_t = end_t - this.trialdur;
        
            idx = find(abs(trs-init_t) <= 1e-1);
            idx_init = idx(1);
            
            % TODO: refine range; maybe throw out the first and last few
            % seconds of data
            range = idx_init + (1:(this.trialdur*this.sfrq));
        
            tt = trs(range);
            xe = xe_rs(range);
            ye = ye_rs(range);

            % Detrend
            if this.condition == 2
                xe = detrend(xe);
                ye = detrend(ye);
            end
            % Align starting point
            this.X = [xe-xe(1); ye-ye(1)]';
            this.t = tt - tt(1);
        end

    end
end


