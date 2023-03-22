classdef handCrossTrialAnalysis < handle
    % Filename:	handCrossTrialAnalysis
    % Authors:  Rika Sugimoto Dimitrova, Federico Tessari, James Hermus
    % Date:		19 Dec 2022
    % Description:	Object oriented program which calls the 
    % handTrialAnalysis object program to computed pooled quantities
    % across all 10 trials of the hand-posture experimenthandTrialAnalysiss.
    
    properties

        subjDex % indices of subjects
        condDex % indices of conditions (1: vision, 2: no vision)
        trialDex % indices of trials

        test % all data
        ensemble % Cell array of structs for each experimental condition
                 % storing hand positions
                 % x: x-coordinate (m) [Subject x trial x N]
                 % y: y-coordinate (m) [Subject x trial x N]
        time % time stamp of data (s) [N x 1]
        
    end
    
    methods
        function this = handCrossTrialAnalysis(subjDex, condDex, trialDex)
            if(nargin==0)
                this.subjDex = 1:10;
                this.condDex = 1:1;
                this.trialDex = 1:10;
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
            for subj = this.subjDex
                for condition = this.condDex
                    for trial = this.trialDex
                        this.test{subj,condition,trial} = ...
                            handTrialAnalysis(subj, condition, trial);
                    end
                end
                disp([int2str(subj),'/',int2str(this.subjDex(end))]);
            end
                        
        end
        
        function [] = get_ensemble(this)
                        
           for condition = this.condDex
               this.time{condition} = ...
                   this.test{this.subjDex(1),condition,this.trialDex(1)}.t;
               isubj = 1;
               for subj = this.subjDex
                   for trial = this.trialDex
                       this.ensemble{condition}.x(isubj,trial,:) = ...
                           unwrap(this.test{subj,condition,trial}.X(:,1));
                       this.ensemble{condition}.y(isubj,trial,:) = ...
                           unwrap(this.test{subj,condition,trial}.X(:,2));
                   end
                   isubj = isubj + 1;
               end
           end

        end
        

    end
end

