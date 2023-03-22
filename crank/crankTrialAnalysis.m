classdef crankTrialAnalysis < handle 
    % Data intial processing for crank-turning data Analysis
    % Filename:	crankTrialAnalysis.m
    % Author:  James Hermus
    % Date:     18 Nov 2022
    % Description: When initalized, this class imports crank-turning
    % subject data and sorts it.
    
    properties
        
        % Subject Specific Parameters
        subject % subject number in study
        trial % Trial Number (Feedback only) % ADD FIX THIS
        speed % Speed and direction are incorporated , according to speedNames
        names % subject identifier
        speedNames % Name of speed + direction conditions
        l1a % Length of upper arm (m)
        l2a % Length of forearm (from elbow to center of fist) (m)
        ch  % Length from wrist to center of fist (m)
        wt  % Subject weight (kg)
        
        % Trial Spesific Parameters
        lc    % Crank Radius (m)
        d1    % Horiz distance from crank center to shoulder center (m) (right when facing robot is positive)
        d2    % Vertical distance from crank center to shoulder center (m) (postive is toward the robot when facing it)
        
        % Assumed Human dynamics
        turnSpeed
        
        sfrq % Sampling frequency
        
        % Measures
        t % Time (s)
        X % Handle position: [Nx2] (m) in cartesian crank coordinates
        thcp  % Crank angular postion (rad)
        thcv  % Crank angular velocity (rad/s)
        thca  % Crank angular acceleration (rad/s^2)
        F_rot % Force normal and tangental [Nx2] (m)
        q     % shoulder and elbow realitive joint angles [Nx2] (rad)
        q_dot % shoulder and elbow realitive joint velocities [Nx2] (rad)

    end
    
    methods
        function  [this] = crankTrialAnalysis(subject, speed, trial, Z_gain_input)
            this.subject = subject;
            this.speed = speed;
            this.trial = trial;
            
            % During construction always run
            this.sfrq = 200; % Enforced sampling frequency
            this.get_subjectParam(); % Defines subject paramters
            this.get_subjectData();
            
        end
        
        %% Main analysis functions:
        function get_subjectParam(this)
            
            % General Shared Variables
            this.names = {'asn','bgn','iwd','onz','iuz','ktn','opo','inz','qba','idg'};
            this.speedNames = {'CW, Slow', 'CW, Med', 'CW, Fast', 'CCW, Slow', 'CCW, Med', 'CCW, Fast'};
            
            turnSpeedArray = [0.075,0.5,2,0.075,0.5,2];
            this.turnSpeed = turnSpeedArray(this.speed);
            
            %% HUMAN:
            [param, grpordr] = this.getparam();
            
            this.l1a = param(1);         % Length of upper arm
            this.l2a = param(2);         % Length of forearm (from elbow to center of fist)
            this.ch = param(3);          % From wrist to center of fist
            this.wt = param(4);          % Weight in kg
            this.lc = param(5);          % length of crank
            this.d1 = param(6);          % Horiz distance from shoulder to crank center
            this.d2 = param(7);          % Vert distance from shoulder to crank center
                        
            
        end
        
        function [] = get_subjectData(this)
            
            [param, grpordr] = this.getparam();
            jdex = find(this.speed == grpordr);
            
            % Select trial number skipping catch trials
            % From Joe's code the blind matrix defines which trials did not
            % have visual feedback and which did. Assumes same pattern of
            % feedback for same speed and direction.
            blind = [0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0;...
                0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0;...
                0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0;...
                0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0;...
                0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0;...
                0,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0];
            
            dex_Trial_Vision = find(~blind(this.speed,:));
            
            % Load data using:
            % jdex to switch to direction speed/direction
            % dex_Trial_Vision to choose only trials with had visual feedback
            % one is added to this.trial to skip the first trial in all cases.
            data = load(['BrownData/crank_data/',...
                this.names{this.subject},'/', this.names{this.subject}, int2str(jdex), sprintf('%02d',dex_Trial_Vision(this.trial + 1)) ,'.ASC']);
            
            % Compute position as positive angle
            p = atan2(sin(data(:,1)), cos(data(:,1)));
            indxs = find(p(:,1)<0);
            p(indxs,1) = p(indxs) + 2*pi;
            
            % Define Fitler Motion
            cf = 10; % Cutoff frequency
            [b,a] = butter(2,cf/(this.sfrq/2));
            
            % Define Data position, velocity, and acceleration
            thcp = p; % (ERROR FIXED - adding pi was an error corrected 3/28/18 ) + pi; % Crank Angular Position
            thcv = filtfilt(b,a, data(:,2)*2*pi ); % Crank Angular Velocity -- FILTERING VELOCITY
            % thca = vel2acc(thcv); % Crank Angular Acceleration
            % thca = diff(thcv)*200; % Crank Angular Acceleration
            % thca = [thca; thca(end)];
            
            %Check filter Motion
            % figure;
            % subplot(2,1,1); hold on; plot(data(:,2)*2*pi); plot(thcv); hold off; title('velocity');
            % subplot(2,1,2); hold on; plot(thca); hold off; title('acceleration');
            
            % Define Filter Normal Forces
            [b,a] = butter(2,cf/(this.sfrq/2));
            F_raw_Norm = data(:,3);
            F_raw_Norm_filt = filtfilt(b,a,F_raw_Norm);
            
            % Define Filter Tangential Forces
            % Tang_cf = [0.5,2,10,0.5,2,10];
            Tang_cf = [0.5,10,10,0.5,10,10];
            cf = Tang_cf(this.speed); % Cutoff frequency
            [b,a] = butter(2,cf/(this.sfrq/2));
            F_raw_Tang = data(:,4);
            F_raw_Tang_filt = filtfilt(b,a,F_raw_Tang);
            
            % Check filter
            % figure;
            % subplot(2,1,1); hold on; plot(F_raw_Norm); plot(F_raw_Norm_filt,'linewidth',1.5); hold off; title('Normal');xlim([1 length(F_raw_Norm_filt)]);
            % xlabel('Sample Number'); ylabel('Force (N)'); set(gca, 'fontsize',14);legend('Raw Signal','Filtered Signal');
            % subplot(2,1,2); hold on; plot(F_raw_Tang); plot(F_raw_Tang_filt,'linewidth',1.5); hold off; title('Tangential');
            % xlabel('Sample Number'); ylabel('Force (N)'); set(gca, 'fontsize',14); legend('Raw Signal','Filtered Signal');xlim([1 length(F_raw_Norm_filt)]);
            
            [ xe, ye, q1p, q2p ] = inverseKinematics(this,thcp);
            t = [0:1/this.sfrq:1/this.sfrq*length(thcp)-1/this.sfrq]';
            
            %% Cut all data to only full cycles
            dexReset = find(abs(diff(thcp)) > 4);
            dexRange = [dexReset(1):dexReset(end)];
            
            this.thcp = thcp(dexRange);
            this.thcv = thcv(dexRange);
            this.F_rot = [F_raw_Norm_filt(dexRange)';F_raw_Tang_filt(dexRange)']; % Normal and Tangential Forces
            this.t = t(dexRange) - t(dexRange(1));
            this.X = [xe(dexRange), ye(dexRange)];
            this.q = [q1p(dexRange), q2p(dexRange)];
            %figure(4*this.subject-3);hold on;plot(t,unwrap(p)-unwrap(p(1)));
            %figure(4*this.subject-2);hold on;plot(this.t,unwrap(this.thcp));
            %figure(4*this.subject-1);hold on;plot(t,unwrap(p)-unwrap(p(1))+this.turnSpeed*2*pi*t);
            %figure(4*this.subject);hold on;plot(this.t,unwrap(this.thcp)+this.turnSpeed*2*pi*this.t);

            
            %             figure;
            %             subplot(3,1,1); plot(this.t, this.thcp); ylabel('Postion (rad)');
            %             subplot(3,1,2); plot(this.t, this.thcv); ylabel('Velocity (rev/s)');
            %             subplot(3,1,3); plot(this.t, this.F_rot); ylabel('Force (N)');
            %             xlabel('Time (sec)');
            
            
        end

        function [ xe, ye, q1p, q2p, thea ] = inverseKinematics( this, thc )
            %Compute joint angles
            %   Computes crank postion and joint angles input params provides the
            %   subject data from the get"kino number"param()
            
            % Compute crank postion in carteasion
            xe = this.d1 + this.lc*cos(thc);
            ye = this.d2 + this.lc*sin(thc);
            ld = sqrt(xe.^2 + ye.^2);
            
            % Right Hand
            alpha1 = atan2(ye,xe);
            alpha2 = acos((xe.^2 + ye.^2 + this.l1a^2 - this.l2a^2)./(2*this.l1a*sqrt(xe.^2 + ye.^2)));
            alpha3 = acos((xe.^2 + ye.^2 + this.l2a^2 - this.l1a^2)./(2*this.l2a*sqrt(xe.^2 + ye.^2)));
            
            % spatial angles
            q1p = alpha1 - alpha2;
            q2p = alpha1 + alpha3 - q1p;
            
            %% Old Code that works
            % % Compute elbow angle relative
            % q2p = pi - acos((-ld.^2+l1a^2+l2a^2)./(2*l1a*l2a));
            %
            % % shoulder angle
            % q1p =  atan2(ye,xe) - acos((ld.^2+l1a^2-l2a^2)./(2.*ld.*l1a));
            
            % elbow angle absolute
            thea = q1p + q2p;
            
        end
        
        function [J] = get_jacobian(this,q)
            
            q1p = q(1);
            q2p = q(2);
            
            J11 = -( this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p) );
            J12 = -this.l2a*sin(q1p + q2p);
            J21 = this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p);
            J22 = this.l2a*cos(q1p + q2p);
            J = [J11, J12; J21, J22];
            
        end
                
        function [ X_bin, vel_bin ] = customBinning(this, x, y, vel, p)
            
            N = 200; % number of bins
            posEdges = linspace(0,2*pi,N+1);
            
            for j = 1:N
                
                binDex = find(p > posEdges(j) & p <= posEdges(j+1)); % Find index pos prosition with in each bin
                inda(j,1:length(binDex)) = binDex;                   % save index ast inda
                numPerBin(j) = length(binDex);
                vel_bin(j) = mean(vel(binDex));                  % take mean of velocity positon points with in each bin
                vel_bin_tot{j} = vel(binDex);
                
                % Find mean radius too
                x_bin(j) = mean(x(binDex));
                y_bin(j) = mean(y(binDex));
                r(j) = sqrt(mean(x(binDex)).^2+mean(y(binDex)).^2);
                
            end
            X_bin = [x_bin; y_bin];
            
        end
        
        %% Get paramerters for specific subjects
        % Recycle structure from Joe
        function [param, grpordr] = getparam(this)
            
            % function [param, grpordr] = getparam(subj)
            %
            % Another way of getting subject parameters -- like a front end to
            % the different get'subj'params.
            
            if(    this.subject == 1 )
                [param, grpordr] = this.getasnparam();
            elseif(this.subject == 2 )
                [param, grpordr] = this.getbgnparam();
            elseif(this.subject == 3 )
                [param, grpordr] = this.getiwdparam();
            elseif(this.subject == 4 )
                [param, grpordr] = this.getonzparam();
            elseif(this.subject == 5 )
                [param, grpordr] = this.getiuzparam();
            elseif(this.subject == 6 )
                [param, grpordr] = this.getktnparam();
            elseif(this.subject == 7 )
                [param, grpordr] = this.getopoparam();
            elseif(this.subject == 8 )
                [param, grpordr] = this.getinzparam();
            elseif(this.subject == 9 )
                [param, grpordr] = this.getqbaparam();
            elseif(this.subject == 10 )
                [param, grpordr] = this.getidgparam();
            else
                error('Subject get params not selected.');
            end
            
        end
        
        function [param, grpordr] = getasnparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,3,4,6,5,1];
            
            la = 14.5*0.0254;   % Length of upper arm
            lf = 14.75*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3.0*0.0254;   % From wrist to center of fist
            wt = 175*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 18.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getbgnparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [6,4,2,5,1,3]; % Omit set '0'
            
            la = 14.25*0.0254;   % Length of upper arm
            lf = 14.5*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 160*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 18.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getidgparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,4,3,6,1,5]; % Omit set '0'
            
            la = 12.75*0.0254;   % Length of upper arm
            lf = 13.0*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 173*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.125*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getinzparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [1,6,5,3,4,2]; % Omit set '0'
            
            la = 14.25*0.0254;   % Length of upper arm
            lf = 15.0*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3.0*0.0254;   % From wrist to center of fist
            wt = 205*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.75*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getiuzparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,4,3,6,5,1]; % Omit set '0'
            
            la = 13.75*0.0254;   % Length of upper arm
            lf = 14*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 160*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getiwdparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [1,3,5,4,6,2]; % Omit set '0'
            
            la = 11.5*0.0254;   % Length of upper arm
            lf = 13.25*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3*0.0254;   % From wrist to center of fist
            wt = 185*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 16.75*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getktnparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [5,1,6,2,3,4]; % Omit set '0'
            
            la = 14.75*0.0254;   % Length of upper arm
            lf = 14.75*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 3.0*0.0254;   % From wrist to center of fist
            wt = 155*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
            
        end
        
        function [param, grpordr] = getonzparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,3,1,4,5,6]; % Omit set '0'
            
            la = 13.5*0.0254;   % Length of upper arm
            lf = 14.5*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 155*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
        function [param, grpordr] = getopoparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [2,1,6,3,4,5]; % Omit set '0'
            
            la = 14.5*0.0254;   % Length of upper arm
            lf = 15.25*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 160*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
            
        end
        
        function [param, grpordr] = getqbaparam(this)
            
            % This function gets the parameters measured off the subject.
            
            grpordr = [1,6,5,3,2,4]; % Omit set '0'
            
            la = 13.5*0.0254;   % Length of upper arm
            lf = 13.25*0.0254;     % Length of forearm (from elbow to center of fist)
            ch = 2.75*0.0254;   % From wrist to center of fist
            wt = 180*0.454;     % Weight in kg
            
            lc = (2.3 + 3.5/2)*0.0254;  % Crank Radius
            d1 = 0.0;                   % Horiz distance to crank center
            d2 = 17.5*0.0254;                  % Vert distance to crank center
            
            param = [la lf ch wt lc d1 d2];
        end
        
    end
end


