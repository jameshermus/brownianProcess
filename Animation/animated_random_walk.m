 clc, close all
set(0, 'DefaultLineLineWidth', 1);
set(groot,'defaultAxesFontSize',15);
set(groot,'defaultAxesColor',[1 1 1]); % axes background
set(0,'defaultfigurecolor',[1 1 1]); % figure background
set(groot,'defaultAxesBox','on');    % box on
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','none');
set(groot,'defaultTextInterpreter','none');

global KEY_IS_PRESSED

KEY_IS_PRESSED = 0;

disp('Press any key to end program')

flag = true;

while(flag == true)

    for subj = randi(10) % select one subject
        %-- upperLimb --%
        if this.bool_analyzeUpperLimb
            trial = randi(10);
            fullscreen = get(0,'ScreenSize');
            f = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
            ax = axes;
            set(gcf, 'KeyPressFcn', @myKeyPressFcn)
            hold on
            %-- upperLimb: vision --%
            %Colors
            grey = [0.7 0.7 0.7];
            orange = [0.95 0.8 0.6];
            red_mat = [0.8500 0.3250 0.0980];
            h1 = animatedline('Color',grey,'LineWidth',2,'LineStyle','-'); %Other Colors
            h2 = animatedline('Color','k','MaximumNumPoints',100,'LineWidth',2,'LineStyle','-');
            p = plot(ax,squeeze(1e3*(this.xe_v(subj,trial,:))),squeeze(1e3*(this.ye_v(subj,trial,:))),'ro','MarkerSize',7,'LineWidth',2);

            xlim([min(squeeze(1e3*(this.xe_v(subj,trial,:)))) max(squeeze(1e3*(this.xe_v(subj,trial,:))))])
            ylim([min(squeeze(1e3*(this.ye_v(subj,trial,:)))) max(squeeze(1e3*(this.ye_v(subj,trial,:))))])
            title(strcat('Subj = ',num2str(subj),' Trial = ',num2str(trial)))
            xlabel('X direction [mm]');
            ylabel('Y direction [mm]');

            for idx = 1:5:length(squeeze(this.xe_v(subj,trial,:)))
                %             tic
                p.XData = 1e3*(this.xe_v(subj,trial,idx));
                p.YData = 1e3*(this.ye_v(subj,trial,idx));
                addpoints(h1,1e3*(this.xe_v(subj,trial,idx)),1e3*(this.ye_v(subj,trial,idx)));%,'Color',[0 0.4470 0.7410]);
                addpoints(h2,1e3*(this.xe_v(subj,trial,idx)),1e3*(this.ye_v(subj,trial,idx)));%,'Color',[0 0.4470 0.7410]);
                legend(p,strcat('Time = ', num2str(floor(this.t_v(idx)))))
                %             toc
                pause(0.020);
                if KEY_IS_PRESSED
                    flag = false;
                    close all
                    break;
                end
            end
        end
    end

end

function myKeyPressFcn(hObject, event)
    global KEY_IS_PRESSED
    KEY_IS_PRESSED  = 1;
    disp('Key pressed: ending program')
end