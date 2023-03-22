classdef simulateModel < handle 
    % Model simulation
    % Filename:	simulateModel.m
    % Author:  Federico Tessari, Rika Sugimoto Dimitrova
    % Date:     19 March 2023
    % Description: When initialized, this class simulates the descriptive
    % model with forward-path velocity control and intermittent position
    % feedback.

    properties

        body_model
        control_mode
        N_trials

        simOut
        
    end

    methods
        function this = simulateModel(body_model, control_mode, N_trials)
            if nargin==0
                prompt = "What system do you want to analyze? Press 1 for Upper Limb, Press 2 for Lower Limb: ";
                body_model = input(prompt);

                if body_model == 1
                    prompt = "Selected control mode: press 1 for Static Posture, Press 2 for Crank Turning: ";
                    control_mode = input(prompt);
                elseif body_model == 2
                    control_mode = 1;
                else
                    error("Invalid model selection.");
                end

                prompt = "How many trials do you want to run?";
                N_trials = input(prompt);
                if N_trials < 1
                    N_trials = 1;
                else
                    N_trials = round(N_trials);
                end                
            end

            % Set parameters in the Simulink model
            if body_model == 1
                m = 1.5; %[kg]
                b = 30; %[Nms/rad]
                g = 0; %[m/s^2]
                L = 0.35; %[m]
                ct = 0; %Control Type
            elseif body_model == 2
                m = 75; %[kg]
                b = 226; %[Nms/rad]
                g = 9.81; %[m/s^2]
                L = 1; %[m]
                ct = 1; %Control Type
            else
                error("Invalid model selection.");
            end%body_model
            f_pole = (b/m)/(2*pi);

            if control_mode == 1
                cm = 1; % Static Posture
                wref = 0; %[rpm]->[rad/s]
                n_gain = 1;
                Tsim = 240;
            elseif control_mode == 2
                cm = 2; % Crank Turning
                wref = 5*(2*pi)/60; %[rpm]->[rad/s]
                n_gain = 7;
                Tsim = 20*(2*pi)/wref;
            else
                error("Invalid control mode");
            end%control_mode
            
            % Simulate model 
            for i = 1:N_trials
                % Motor noise seed
                mns = randi(1e5,1);
                % Disturbance noise seed
                dns = randi(1e5,1);
                % Sensor noise seed
                sns = randi(1e5,1);
                % Neural noise seed
                nns = randi(1e5,1);

                modelname = 'sim_ff_velocity';
                simIn = Simulink.SimulationInput(modelname);
                simIn = setVariable(simIn,'mns',mns);
                simIn = setVariable(simIn,'dns',dns);
                simIn = setVariable(simIn,'sns',sns);
                simIn = setVariable(simIn,'nns',nns);
                simIn = setVariable(simIn,'wref',wref);
                simIn = setVariable(simIn,'n_gain',n_gain);
                simIn = setVariable(simIn,'cm',cm);
                simIn = setVariable(simIn,'m',m);
                simIn = setVariable(simIn,'b',b);
                simIn = setVariable(simIn,'g',g);
                simIn = setVariable(simIn,'L',L);
                simIn = setVariable(simIn,'ct',ct);
                simIn = setModelParameter(simIn,'StopTime',num2str(Tsim));

                this.simOut{i} = sim(simIn);
            end%i

            this.body_model = body_model;
            this.control_mode = cm;
            this.N_trials = N_trials;

        end%simulateModel


    end%methods

end