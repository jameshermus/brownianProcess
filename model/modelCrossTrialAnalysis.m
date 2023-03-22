classdef modelCrossTrialAnalysis < handle
    % Filename:	modelCrossTrialAnalysis
    % Author:  Rika Sugimoto Dimitrova
    % Date:		19 March 2023
    % Description:	Object oriented program which assembles the ensemble
    % for the simulation data
    
    properties

        sysDex
        trialDex
        datafiles

        time
        ensemble % Cell array of simulation data
        
    end
    
    methods
        function this = modelCrossTrialAnalysis(sysDex, trialDex, files)
            if(nargin==0)
                this.sysDex = 1:3; % 1: crank turning, 2: hand posture, 3: quiet standing
                this.trialDex = 1:200;
                this.datafiles = {'out_v_crank','out_v_hand','out_v_stand'};
            else
                this.sysDex = sysDex;
                this.trialDex = trialDex;
                for i = 1:length(files)
                    this.datafiles{sysDex(i)} = files{i};
                end
            end
                        
            this.get_ensemble();
            
        end
        
        function [] = get_ensemble(this)

            %Loading Simulation Data
            if(~license('test','simulink'))
                error('It was detected that Simulink is not currently installed. In order to run the process_models function please install simulink.');
            end            
                     
            for system = this.sysDex
                clear simData
                simData = load(['BrownData/modelling/' this.datafiles{system} '.mat']);
                for trial = this.trialDex
                    this.ensemble{system}.pos(1,trial,:) =...
                        simData.out{trial}.x3.Data - simData.out{trial}.x3.Data(1);
                end
                this.time{system} = simData.out{trial}.x3.Time;
            end


        end
        

    end
end

