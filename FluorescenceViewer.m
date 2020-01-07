classdef FluorescenceViewer < handle
    
    properties
        data = [];
        bgclr = 'w';
        txtclr = 'k';
        ChCLR
        XRange = [350, 900];
        Main = '';
        FNMS
        Table
        t
        AxesAllFluor
        Axes
    end
    methods
        % Build the ATLAS graphical user interface
        function app = FluorescenceViewer
            CurrentPath = fileparts(mfilename('fullpath'));
            app.Main = [CurrentPath filesep 'Spectra Raw Data'];
            if ~isfolder(app.Main)
                app.Main = uigetdir(pwd, 'Select path to fluoresence data');
            end
            addpath(CurrentPath);
            app.FNMS = dir('*.csv*');
            app.FNMS = struct2table(app.FNMS);
            app.FNMS = app.FNMS.name;
            app.FNMS = strrep(app.FNMS, '.csv', '');
            app.FNMS = [{' '}, app.FNMS(:)'];
            func_newfig(app);
            fig = gcf;
            func_plot(app, [ 'fig' num2str(fig.Number)]);
            
        end % end of Main function
        
        
        % Make a new figure (Includes gating buttons)
        function func_newfig(app)
            %% Open a new figure for the data table
            fig = figure;
            fig.Color = app.bgclr;
            fig.InvertHardcopy = 'off';
            fig.Position = [10, 10, 640 500];
            %% Put in a table
            
            % Put in some default values
            DefaultABS = {'SIRPa', 'CD11c', 'CD8', 'MHC-II', 'CLEC9a', 'PCREB', 'IRF4', 'B220','CD64' 'CD3', 'PD1', 'Foxp3', 'PD-L1'};
            DefaultFLR = {'BV421', 'BV480', 'BV510', 'DY395XL', 'CF633', 'AF488', 'AF514', 'Atto490LS', 'CF660C', 'AF700', 'CF750', 'eFluor570', 'CF594'};
            DefaultDET = {'418-446', '460-486', '501-548', '583-620', '649-673', '498-521', '536-562', '586-697', '672-695', '704-745', '761-800', '569-592', '615-647'};
            DefaultLAS = {405,405,405,405,633,489,489,489,662,662,662,561,600};
            DefaultDETType = {'PMT1','HyD2','PMT3','HyD4','PMT5','HyD2','HyD4','PMT5','HyD2','HyD4','PMT5','HyD2','HyD4'};

            DefaultCHN = {1};
            DefaultSEQ = {1,1,1,1,1,2,2,2,3,3,3,4,4};
            
            app.Table = cell(13,9);
            app.Table(1:13,1)= {true};
            
            app.Table(1:13,2)= DefaultCHN;
            app.Table(1:13,3)= DefaultABS;
            app.Table(1:13,4)= DefaultFLR;
            app.Table(1:13,5)= DefaultSEQ;
            
            app.Table(1:13,6)= DefaultDET;
            app.Table(1:13,7)= DefaultLAS;
            app.Table(1:13,8)= {1};
            app.Table(1:13,9)= DefaultDETType;

            
            % Create the table of options
            app.t = uitable(fig);
            
            % This makes drop down menus right in the table
            app.t.ColumnFormat = ({[] [] [] app.FNMS [] [] [] []});
            app.t.Data = app.Table;
            app.t.Position = [0 0 1500 500];
            app.t.ColumnName = {'Plot', 'Weight', 'Antibody', 'Fluorophore', 'Sequential', 'Detector', 'Laser', 'Laser Power', 'Detector Type'};
            app.t.ColumnEditable = true;
            app.t.ColumnWidth = {50, 50, 70, 70, 70, 52, 50, 75, 80};
            
            uicontrol('Style', 'pushbutton', 'String', 'Save',...
                'Position', [10 10 50 30],...
                'Callback', @(~, ~) func_save(app));
            
            uicontrol('Style', 'pushbutton', 'String', 'Load',...
                'Position', [70 10 50 30],...
                'Callback', @(~, ~) func_load(app));
            
            uicontrol('Style', 'pushbutton', 'String', 'Load Comp Matrix',...
                'Position', [140 10 150 30],...
                'Callback', @(~, ~) func_CmpMatLoad(app));
            
            uicontrol('Style', 'pushbutton', 'String', 'Plot Comp Matrix',...
                'Position', [300 10 150 30],...
                'Callback', @(~, ~) func_CmpMatPlot(app));
            
            %% Open a new figure for the Spectra
            fig = figure;
            fig.Color = app.bgclr;
            fig.InvertHardcopy = 'off';
            fig.Position = [10, 10, 1000 600];
            app.data.FigNum = [ 'fig' num2str(fig.Number)]; %(Use this to assign handles)
            
            app.t.CellEditCallback = @(~,event) func_plot(app, app.data.FigNum);

        end     % end of new figure function 
        
        % Plot Data function
        function func_plot(app, num)
            %% re-order the list of fluorophores based on ch number
            
            logic = [app.t.Data{:, 1}]==1;
            subdat = app.t.Data(logic, :);
            
            [~, ind] = sort([subdat{:, 5}]);
            % Pull the detector ranges of the sorted sequentials
            RNG = subdat(ind, 6);
            ind2 = cell2mat(strfind(RNG, '-'));
            for i=1:numel(RNG)
                RNG{i} = RNG{i}(1:ind2(i)-1);
            end
            RNG = str2double(RNG);
            % Sort each equential
            for n=1:numel(unique([subdat{:, 5}]))
                SubRNG = RNG([subdat{ind, 5}]==n);
                % Find the sorted index for this sequential
                [~, Subind] = sort(SubRNG);
                SubI2 = ind([subdat{ind, 5}]==n); 
                ind([subdat{ind, 5}]==n) = SubI2(Subind);
            end
            subdat = subdat(ind, :);

            fig = gcf;
            % Select the figure with the spectral plots
            figure(str2double(strrep(num, 'fig', '')));
            clf
            %% Initialize plots
            
            % Pull the number of sequentials and channels from the table
            NSeq = max([subdat{:, 5}]);
            Fluors = subdat(:, 4);
            Nch = numel(Fluors);
            ChMat = zeros(numel([subdat{:, 5}]), numel([subdat{:, 5}]));
            % Calculate subplot dimensions
            SubPltD1 = ceil(sqrt(NSeq))+2;
            SubPltD2 = ceil(sqrt(NSeq));
            
            app.AxesAllFluor.(num) = subplot(SubPltD1, SubPltD2, 1:SubPltD2);
            app.AxesAllFluor.(num).XLabel.String = 'Wavelength, nm';
            app.AxesAllFluor.(num).XLim = app.XRange;
            app.AxesAllFluor.(num).YLim = [0 1.2];
            app.AxesAllFluor.(num).YColor = app.bgclr;
            app.AxesAllFluor.(num).Color = app.bgclr;
            app.AxesAllFluor.(num).XColor = app.txtclr;
            app.AxesAllFluor.(num).YLabel.String = 'All Spectra';
            app.AxesAllFluor.(num).YLabel.Color = app.txtclr;
            hold on
            %% Initialize other plots
            for i = 1:NSeq
                SubFigN = ['N' num2str(i)];
                app.Axes.(num).(SubFigN) = subplot(SubPltD1, SubPltD2, SubPltD1*SubPltD2-i-1);
                app.Axes.(num).(SubFigN).XLabel.String = 'Wavelength, nm';
                app.Axes.(num).(SubFigN).XLim = app.XRange;
                app.Axes.(num).(SubFigN).YLim = [0 1.25];
                app.Axes.(num).(SubFigN).YColor = app.bgclr;
                app.Axes.(num).(SubFigN).Color = app.bgclr;
                app.Axes.(num).(SubFigN).XColor = app.txtclr;
                app.Axes.(num).(SubFigN).YLabel.String = ['Sequential ' num2str(NSeq-i+1)];
                app.Axes.(num).(SubFigN).YLabel.Color = app.txtclr;
                hold on
            end
            %% Plot on the currently selected figure         
            
            % Clear the current axes
            cla(app.AxesAllFluor.(num))
            for i = 1:NSeq
                SubFigN = ['N' num2str(i)];
                cla(app.Axes.(num).(SubFigN))
            end
            
            % Build the color pallet
            clr = jet;
            ind = 1:round(size(clr, 1)/(numel([subdat{:, 2}]))):size(clr, 1);
            if numel(ind)~=numel([subdat{:, 2}])
                ind = 1:round(size(clr, 1)/(numel([subdat{:, 2}])+2)):size(clr, 1);
                ind = ind(1:numel([subdat{:, 2}]));
            end
            clr = clr(ind, :);
            Allclr = clr;
            clr = clr([subdat{:, 1}], :);
            
            % Load Relevent Spectra
            cd(app.Main);
            Weights = [subdat{[subdat{:, 1}], 2}];
            LaserWeights = [subdat{[subdat{:, 1}], 8}];
            for n=1:Nch
                FileNames = dir(['*' Fluors{n} '.csv*']);
                FileNames=struct2table(FileNames);
                FileNames=FileNames.name;
                dat=importdata(FileNames);
                datDYE{n} = dat.data;
                % Normalize the data
                datDYE{n}(:,2) = datDYE{n}(:,2)/max(datDYE{n}(:,2));
                datDYE{n}(:,4) = datDYE{n}(:,4)/max(datDYE{n}(:,4));
                % Multiply by Weights
                datDYE{n}(:,4) = datDYE{n}(:,4).*Weights(n);
            end
            
            % Plot the spectra
            for n=1:Nch
                % Treat spectra and laser power differently
                if sum(strcmp(Fluors{n}, {'2PPower_Chameleon','PMT_Efficiency','HyD_Efficiency'}))~=0
                    % Plot the  spectra   
                    plot(app.AxesAllFluor.(num),datDYE{n}(:,1),datDYE{n}(:,2) ...
                        ,'--', 'Color', app.txtclr, 'Linewidth', 0.5, 'DisplayName',Fluors{n})
                    % Label the plot
                    [~, ind] = max(datDYE{n}(:,4));
                    txt = text(app.AxesAllFluor.(num), datDYE{n}(ind,3), 1.18, Fluors{n});
                    txt.Color = app.txtclr;
                else
                    % Plot the Emission spectra
                    plot(app.AxesAllFluor.(num),datDYE{n}(:,3),datDYE{n}(:,4)/max(datDYE{n}(:,4)), ...
                        'Color', clr(n,:), 'Linewidth', 1.5, 'DisplayName',[Fluors{n} ', EM'])            
                    % Plot the Absorption spectra   
                    plot(app.AxesAllFluor.(num),datDYE{n}(:,1),datDYE{n}(:,2) ...
                        ,'--', 'Color', clr(n,:), 'Linewidth', 0.5, 'DisplayName',[Fluors{n} ', EX'])
                    % Label the plot
                    [~, ind] = max(datDYE{n}(:,4));
                    txt = text(app.AxesAllFluor.(num), datDYE{n}(ind,3), 1.18, [Fluors{n} ':' app.t.Data{n, 3}]);
                    txt.Color = clr(n,:);
                end
                txt.FontSize = 7;
                txt.HorizontalAlignment = 'center';
                txt.Rotation = 45;
            end % end plotting all spectra
            
            % Find the number of sequentials
            NSeq = max([subdat{:, 5}]);
            
            % Plot the spectra in their respective sequential
            for n=1:Nch
                % Treat spectra and laser power differently
                if sum(strcmp(Fluors{n}, {'2PPower_Chameleon','PMT_Efficiency','HyD_Efficiency'}))~=0
                    continue
                end
                    
                % Pull out the sequential #
                i = [subdat{:, 5}];
                i = i(n);
                
                % Pull out laser Wavelength
                las = [subdat{:, 7}];
                las = las(n);
                % Find the index of the nearest wavelength
                [~, ind] = min(abs(datDYE{n}(:,1)-las));
                relEM = datDYE{n}(ind,2).*LaserWeights(n);
                
                SubFigN = ['N' num2str(NSeq-i+1)];
                % Plot the absorption spectra
                plot(app.Axes.(num).(SubFigN),datDYE{n}(:,3),datDYE{n}(:,4).*relEM, ...
                    'Color', clr(n,:), 'Linewidth', 1.5, 'DisplayName',[Fluors{n} ', EM'])            
                % Plot the excitation spectra   
                plot(app.Axes.(num).(SubFigN),datDYE{n}(:,1),datDYE{n}(:,2) ...
                    ,'--', 'Color', clr(n,:), 'Linewidth', 0.5, 'DisplayName',[Fluors{n} ', EX'])
            end % end plotting spectra in Sequentials
           
            % Plot the laser and detectors
            for n=1:numel([subdat{:, 5}])
                % Pull out the sequential #
                i = [subdat{:, 5}];
                i = i(n);
                SubFigN = ['N' num2str(NSeq-i+1)];
                if sum(strcmp([subdat{n, 4}], {'2PPower_Chameleon','PMT_Efficiency','HyD_Efficiency'}))==0
                    % Pull the detector ranges from the table
                    RNG = subdat(:, 6);
                    RNG = RNG{n};
                    % Convert from text to double
                    st = str2double(RNG(1:(strfind(RNG, '-')-1)));
                    en = str2double(RNG((strfind(RNG, '-')+1):end));
                    % Plot the detector ranges
                    plot(app.Axes.(num).(SubFigN),[st, st, en, en],[0, 1.1, 1.1, 0], ...
                        ':', 'Color',app.txtclr, 'LineWidth', 0.5)
                    txt = text(app.Axes.(num).(SubFigN),(st+((en-st)/2)), 1.18, ['    ' num2str(n) ':' subdat{n, 3}]);
                    txt.FontSize = 7;
                    txt.HorizontalAlignment = 'center';
                    txt.Rotation = 45;
                    txt.Color = app.txtclr;

                    % Pull the Laser excitation from the table
                    LAS = [subdat{:, 7}];
                    LAS = LAS(n);
                    % Plot the detector ranges
                    plot(app.Axes.(num).(SubFigN),[LAS, LAS],[0, 1.1], ...
                        app.txtclr, 'LineWidth', 1)
                end
            end % end plotting laser and detectors
            
            % Plot the scaled parasitic spectra in each sequential
            for i=1:NSeq
                SubFigN = ['N' num2str(NSeq-i+1)];
                % Pull out the laser wavelengths
                [las, ia] = unique([subdat{[subdat{:, 5}]==i, 7}]);
                LsW = LaserWeights([subdat{:, 5}]==i);
                LsW = LsW(ia);
                for n=1:Nch
                    % Treat spectra and laser power differently
                    if sum(strcmp(Fluors{n}, {'2PPower_Chameleon','PMT_Efficiency','HyD_Efficiency'}))~=0
                        % Plot the  spectra   
                        plot(app.Axes.(num).(SubFigN),datDYE{n}(:,1),datDYE{n}(:,2) ...
                            ,'--', 'Color', app.txtclr, 'Linewidth', 0.5, 'DisplayName',Fluors{n})
                        % Label the plot
                        [~, ind] = max(datDYE{n}(:,4));
                        txt = text(app.Axes.(num).(SubFigN), datDYE{n}(ind,3), 1.18, Fluors{n});
                        txt.Color = app.txtclr;
                        txt.FontSize = 7;
                        txt.HorizontalAlignment = 'center';
                        txt.Rotation = 45;
                        continue
                    end
                    % Pull out the sequential #
                    j = [subdat{([subdat{:, 1}]), 5}];
                    j = j(n);
                    for k=1:numel(las)
                        %Find the index of the nearest wavelength
                        [~, ind] = min(abs(datDYE{n}(:,1)-las(k)));
                        if ind>0
                            relEM = datDYE{n}(ind,2).*LsW(k);
                            if isnan(relEM)
                                relEM = 0;                                      
                            end
                            if relEM > 0.01 && j ~= i
                                plot(app.Axes.(num).(SubFigN),las,relEM,'xm', 'MarkerSize', 5)
                                plot(app.Axes.(num).(SubFigN),datDYE{n}(:,1),datDYE{n}(:,2), ':', 'Color', [0.8,0.8,0.8],  'Linewidth', 0.5)
                            end
                        else
                            relEM = 0;
                        end
                        % Sum up the bleed through
                        if k==1
                            Bleed = datDYE{n}(:,4).*relEM;
                        else
                            Bleed = Bleed + datDYE{n}(:,4).*relEM;
                        end
                        % remove the NaN values
                        Bleed(isnan(Bleed))=0;
                        if k==numel(las) && j ~= i
                            plot(app.Axes.(num).(SubFigN),datDYE{n}(:,3),Bleed,'m', 'Linewidth', 0.5)
                        end % end if
                    end % end for 1 to number of lasers in sequential
                    %% Find the integrated signal in each detector
                    
                    % index of chanels in this sequential
                    ind = [subdat{:, 5}]==i;
                    RNG = subdat(ind, 6);
                    Chls = find([subdat{:, 1}]);
                    for k=1:sum(ind)
                        % Pull the detector ranges from the table
                        SubRNG = RNG{k};
                        % Convert from text to double
                        st = str2double(SubRNG(1:(strfind(SubRNG, '-')-1)));
                        en = str2double(SubRNG((strfind(SubRNG, '-')+1):end));
                        %Find the index of the nearest wavelength
                        [~, indst] = min(abs(datDYE{n}(:,3)-st));
                        %Find the index of the nearest wavelength
                        [~, inden] = min(abs(datDYE{n}(:,3)-en));
                        
                        % integrate fluorescence over the detector range
%                         IntFluor = sum(Bleed(indst:inden));
                        % Fit with a Spline
                        Wavelengths = datDYE{n}(indst:inden,3);
                        Fluorescence = Bleed(indst:inden);
                        if numel(Wavelengths)==1
                            IntFluor = Fluorescence;
                        else
                            for j = 1:numel(Wavelengths)
                                if j==1
                                    IntFluor = Fluorescence(j);
                                else
                                    IntFluor = IntFluor + Fluorescence(j)*(Wavelengths(j-1)-Wavelengths(j));
                                end
                            end
                        end
% % %                         % Find coefficients for spline interpolant
% % %                         spl = spline(datDYE{n}(:,3),datDYE{n}(:,4));
% % %                         splint = fnint(spl);
% % %                         IntFluor = diff(fnval(splint,[st,en]))

                        
                        % put that in the corresponding location in the
                        % overlap matrix ChMat = zeros(Nch, Nch);
                        ChMat(Chls(n),(find(ind, 1)+k-1)) = IntFluor;
                    end % end integrating fluorescence in each detector
                end % end run through each channel
            end % end plot scaled parasitic spectra
            %% Initialize Box plot
            app.Axes.(num).BoxPlot = subplot(SubPltD1, SubPltD2, (SubPltD1*SubPltD2-SubPltD2+1):(SubPltD1*SubPltD2));

            % normalize the spillover matrix
            ChMat = ChMat./sum(ChMat);
            app.data.ChMat = ChMat';
            app.ChCLR = Allclr;
            
            plt = bar(app.Axes.(num).BoxPlot,ChMat',...
                'stacked', 'EdgeColor', 'black', 'FaceColor','flat');
            % Put in the channel colors
            for n=1:size(ChMat, 1)
                plt(n).CData = Allclr(n,:);
            end

            box off
            app.Axes.(num).BoxPlot.XLabel.String = 'Channel';
            app.Axes.(num).BoxPlot.XLim = [0 size(ChMat, 1)+1];
            app.Axes.(num).BoxPlot.YLim = [0 1];
            app.Axes.(num).BoxPlot.Color = app.bgclr;
            app.Axes.(num).BoxPlot.YColor = app.txtclr;
            app.Axes.(num).BoxPlot.XColor = app.txtclr;
            app.Axes.(num).BoxPlot.XTickLabels = subdat(:, 3);
            app.Axes.(num).BoxPlot.XTickLabelRotation = 90;
            app.Axes.(num).BoxPlot.YLabel.String = 'Channel Composition';
            app.Axes.(num).BoxPlot.YLabel.Color = app.txtclr;
           

            %% return control back to the table
            figure(fig);
        end % end of plot function
        
        function func_save(app)
            dat = cell2table(app.t.Data);
            dat.Properties.VariableNames = {'Plot', 'Weight', 'Antibody', 'Fluorophore', 'Sequential', 'Detector', 'Laser', 'LaserPower','DType'};
             % Get the file path
            [file,path] = uiputfile('*.csv','Save file name');
            cd(path);
            writetable(dat,file, 'WriteVariableNames', true);
        end
        
        function func_load(app)
            % Get the file
            [file,path] = uigetfile('*.csv','Select File');
            dat = readtable([path file], 'EmptyValue', 0);
            
            dat = table2cell(dat);
            dat([dat{:, 2}]'==0, :) = [];
            dat(:, 1) = {true};
            app.t.Data = dat;
            func_plot(app, app.data.FigNum);
        end
        
        % This Loads the Comp Matrix
        function func_CmpMatLoad(app)
            [file,path] = uigetfile('*.sdm','Select File');
            cd(path)
            fileID = fopen([path file]);
            app.data.CmpMat = textscan(fileID, '%s', 'Delimiter', '<');
            fclose(fileID);
            app.data.CmpMat = app.data.CmpMat{1};
            app.data.CmpMat = app.data.CmpMat(5:end);
            % dat = dat(:, 4:end);
            app.data.CmpMat = app.data.CmpMat(~contains(app.data.CmpMat, {'/Ch'}));          
            app.data.CmpMat = app.data.CmpMat(~contains(app.data.CmpMat, {'Dye'}));
            IND = strfind(app.data.CmpMat, '>');
            IND = cell2mat(IND);
            app.data.NCh = str2double(app.data.CmpMat{end}(3:(IND(end)-1)));
            for i=1:numel(app.data.CmpMat)
                app.data.CmpMat{i} = str2double(app.data.CmpMat{i}((IND(i)+1):end));
            end
            app.data.CmpMat = cell2mat(app.data.CmpMat);
            app.data.CmpMat = reshape(app.data.CmpMat, [app.data.NCh app.data.NCh])';
            % dat = dat./sum(dat, 2);
            
           func_CmpMatPlot(app) 
        end
        
        % This plots the compensation matrix
        function func_CmpMatPlot(app)

            fig = figure(3);
            clf
            fig.Color = app.bgclr;
            fig.InvertHardcopy = 'off';
            %% Plot the measured compensation matrix
            subplot(2,4,1:2)
            plt = bar(app.data.CmpMat./sum(app.data.CmpMat, 2), 'stacked', 'EdgeColor', 'black', 'FaceColor','flat');
            for n=1:app.data.NCh
                plt(n).CData = app.ChCLR(n,:);
            end
            box off
            ax = gca;
            ax.XLim = [0 app.data.NCh+1];
            ax.YLim = [0 1];
            ax.Color = app.bgclr;
            ax.YColor = app.txtclr;
            ax.XColor = app.txtclr;
            ax.XTickLabels = [];
            ax.YLabel.String = 'Measured Channel Composition';
            ax.YLabel.Color = app.txtclr;
            
            % Plot the heatmap of values
            subplot(2,4,3)
            hm = heatmap(app.data.CmpMat);
            hm.Colormap = redbluecmap;
            hm.XDisplayLabels = app.t.Data(:, 3);
            hm.YDisplayLabels = app.t.Data(:, 3); 
            
            
            %% Plot the theoretical compensation matrix
            subplot(2,4,5:6)
            plt = bar(app.data.ChMat, 'stacked', 'EdgeColor', 'black', 'FaceColor','flat');
            for n=1:app.data.NCh
                plt(n).CData = app.ChCLR(n,:);
            end
            box off
            ax = gca;
            ax.XLabel.String = ' Channel';
            ax.XLim = [0 app.data.NCh+1];
            ax.YLim = [0 1];
            ax.Color = app.bgclr;
            ax.YColor = app.txtclr;
            ax.XColor = app.txtclr;
            ax.XTickLabels = app.t.Data(:, 3);
            ax.XTickLabelRotation = 90;
            ax.YLabel.String = 'Theoretical Channel Composition';
            ax.YLabel.Color = app.txtclr;
            
            % Plot the hearmap of values
            subplot(2,4,7)
            hm = heatmap(app.data.ChMat);
            hm.Colormap = redbluecmap;
            hm.XDisplayLabels = app.t.Data(:, 3);
            hm.YDisplayLabels = app.t.Data(:, 3); 
            
            % Plot the error
            subplot(2,4,4)
%             error = (app.data.CmpMat./sum(app.data.CmpMat, 2)-app.data.ChMat)./(app.data.CmpMat./sum(app.data.CmpMat, 2));
            error = (app.data.CmpMat./sum(app.data.CmpMat, 2)-app.data.ChMat);
            hm = heatmap(error);
            hm.Colormap = redbluecmap;
            hm.XDisplayLabels = app.t.Data(:, 3);
            hm.YDisplayLabels = app.t.Data(:, 3);
            title('measured - theoretical')
            subplot(2,4,8)
            hm = heatmap(error./(app.data.CmpMat./sum(app.data.CmpMat, 2)));
            hm.Colormap = redbluecmap;
            hm.ColorLimits = [-10 10];
            hm.XDisplayLabels = app.t.Data(:, 3);
            hm.YDisplayLabels = app.t.Data(:, 3);
            title('(measured - theoretical)/measured')

        end
    end % end of methods
end % end of class definition