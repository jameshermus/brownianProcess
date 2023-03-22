classdef standCrossTrialAnalysis < handle
    % Preprocessing for standing balance data (Santos et. al. 2017)
    % Filename:	balanceCrossTrialAnalysis
    % Authors:  Rika Sugimoto Dimitrova, James Hermus, Federico Tessari
    % Date:		13 Jan 2023
    % Description:	Object oriented program which calls standTrialAnalysis
    % to computed pooled quantities across all subjects and trials
    % for the quiet standing experiments.
    
    properties

        subjDex % indices of subjects
        condDex % indices of conditions 
                % The four experimental conditions are:
                % 1: 'OR', 2: 'CR', 3: 'OF', 4: 'CF' 
                % O = open eyes, C = closed eyes, R = ridged floor, F = foam floor
        trialDex % indices of trials

        test % all data
        ensemble % Cell array of structs for each experimental condition
                 % storing COM position and joint angle data
                 % xcom:    x-coord of COM (m) [Subject x trial x N]
                 % ycom:    y-coord of COM (m) [Subject x trial x N]
                 % zcom:    z-coord of COM (m) [Subject x trial x N]
        time % time stamp of data (s) [N x 1]
        
    end
    
    methods
        function this = standCrossTrialAnalysis(subjDex, condDex, trialDex)
            if(nargin==0)
                this.subjDex = [1:12 14:17 20 21 23 25:27 40 46:48]; % young subjects with no disabilities 
                                      % 11 females; 15 males
                                      % mean height = 1.71 m
                                      % mean mass   = 71.0 kg
                % this.subjDex = [13 18 19 22 24 28:32 35:36 39]; % elderly subjects with no disabilities 
                this.condDex = 1;
                this.trialDex = 1:3;
            else
                this.subjDex = subjDex;
                this.condDex = condDex;
                this.trialDex = trialDex;
            end
                        
            this.get_testStruct();
            this.get_ensemble();
            
        end
        
        function [] = get_testStruct(this)

            clear test
            this.test = cell(length(this.subjDex),length(this.condDex),length(this.trialDex));
            subj_count = 1;
            for subj = this.subjDex
                for condition = this.condDex
                    for trial = this.trialDex
                        this.test{subj_count,condition,trial} = ...
                            standTrialAnalysis(subj, condition, trial);
                    end
                end
                disp([int2str(subj_count),'/',int2str(length(this.subjDex))]);
                subj_count = subj_count + 1;
            end
                        
        end
        
        function [] = get_ensemble(this)
           
           for condition = this.condDex
               this.time{condition} = ...
                   this.test{1,condition,this.trialDex(1)}.t;
               for subj = 1:size(this.test,1)
                   for trial = this.trialDex
                       this.ensemble{condition}.xcom(subj,trial,:) = ...
                           unwrap(this.test{subj,condition,trial}.X(:,1));
                       this.ensemble{condition}.ycom(subj,trial,:) = ...
                           unwrap(this.test{subj,condition,trial}.X(:,2));
                       this.ensemble{condition}.zcom(subj,trial,:) = ...
                           unwrap(this.test{subj,condition,trial}.X(:,3));
                   end
               end
           end

        end
        

    end
end

