classdef crankCrossTrialAnalysis < handle
    % Filename:	crankCrossTrialAnalysis
    % Author:  James Hermus
    % Date:		18 November 2022
    % Description:	Object oriented program which calls the 
    % crankTrialAnalysis object program to computed pooled quantities
    % across all 21 trials of the crank-turning experiments.
    
    properties

        subjDex
        speedDex
        trialDex

        test % all data
        ensemble % Cell array of 
        
    end
    
    methods
        function this = crankCrossTrialAnalysis(subjDex, speedDex, trialDex)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            if(nargin==0)
                this.subjDex = 1:10;
                this.speedDex = 1; % only slow speed CW trial %1:6;
                % Note speedDex can only be 1 or 1:6
                this.trialDex = 1:21;
            else
                this.subjDex = subjDex;
                this.speedDex = speedDex;
                this.trialDex = trialDex;
            end
                        
            this.get_testStruct();
            this.get_ensemble();
            
        end
        
        function [] = get_testStruct(this)

            clear test
            this.test = cell(length(this.subjDex),length(this.speedDex),length(this.trialDex));
            for subj = this.subjDex
                for speed = this.speedDex
                    for trial = this.trialDex
                        this.test{subj,speed,trial} = crankTrialAnalysis(subj, speed, trial);
                    end
                end
                disp([int2str(subj),'/',int2str(this.subjDex(end))]);
            end
                        
        end
        
        function [depMeas] = remove_outliers(this,depMeas)
        
            % Take out major outliers
            % Subject 5 outliers all over the place CCW medium
            % Subject 8 outliers all over the place CCW medium
            
            % Subject 7 outlier CW medium: 3, 6, CCW medium: 14, 21
            dexNoOutlier = find((this.trialDex ~= 3 )&( this.trialDex ~= 6));
            depMeas(3,1,1,7) = mean(squeeze(depMeas([dexNoOutlier],1,1,7)));
            depMeas(6,1,1,7) = mean(squeeze(depMeas([dexNoOutlier],1,1,7)));
            
            dexNoOutlier = find((this.trialDex ~= 14 )&( this.trialDex ~= 21));
            depMeas(14,1,2,7) = mean(squeeze(depMeas([dexNoOutlier],1,2,7)));
            depMeas(21,1,2,7) = mean(squeeze(depMeas([dexNoOutlier],1,2,7)));
            
            % Subject 9 outliers CW 7,8,11, CCW Medium 12
            dexNoOutlier = find((this.trialDex ~= 7 )&( this.trialDex ~= 8)&( this.trialDex ~= 11));
            depMeas(7,1,1,9) = mean(squeeze(depMeas([dexNoOutlier],1,1,9)));
            depMeas(8,1,1,9) = mean(squeeze(depMeas([dexNoOutlier],1,1,9)));
            depMeas(11,1,1,9) = mean(squeeze(depMeas([dexNoOutlier],1,1,9)));
            
            dexNoOutlier = find((this.trialDex ~= 12 ));
            depMeas(12,1,2,9) = mean(squeeze(depMeas([dexNoOutlier],1,2,9)));
            
            % Subject 10 outliers CCW medium 14, CCW fast 6
            dexNoOutlier = find((this.trialDex ~= 14 ));
            depMeas(14,1,2,10) = mean(squeeze(depMeas([dexNoOutlier],1,2,10)));
            
            dexNoOutlier = find((this.trialDex ~= 6 ));
            depMeas(6,2,2,10) = mean(squeeze(depMeas([dexNoOutlier],2,2,10)));
        
        end
        
        function [] = get_ensemble(this)
            
            minLength = zeros(length(this.speedDex),1);
            for speed = this.speedDex
                for subj = this.subjDex
                    for trial = this.trialDex
                        l(subj,speed,trial) = length(this.test{subj,speed,trial}.thcp);
                    end
                end
            end
            
            l_min = min(min(l,[],3));
            clear l
            if length(l_min) == 1
                l = l_min;
                if(mod(l,2)==0)
                    l = l-1; % Make odd number
                end

            else
                l = min([l_min(1:3);l_min(4:6)])-2;
                for speed = 1:3
                    if(mod(l(speed),2)==0)
                        l(speed) = l(speed)-1; % Make odd number
                    end
                end
            end
                      
           dexL = [1,2,3,1,2,3];
                        
           for speed = this.speedDex
               isubj = 1;
               for subj = this.subjDex
                   for trial = this.trialDex
                       this.ensemble{speed}.thcp_r(isubj,trial,:) = ...
                           unwrap(this.test{subj,speed,trial}.thcp(2:l(dexL(speed))));
                   end
                   isubj = isubj + 1;
               end
           end

        end
        

    end
end

