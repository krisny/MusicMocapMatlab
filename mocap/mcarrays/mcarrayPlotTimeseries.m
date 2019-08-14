function mcarrayPlotTimeseries(varargin)
% 
% Plotting all entries in the array. 
% Only suitable for MoCap Data or Norm Data types.
% Based on mcplottimeseries from Mocap Toolbox v1.5
% Adapted by Kristian Nymoen, RITMO/University of Oslo, 2019
%
% 
% syntax
% mcarrayPlotTimeseries(d, marker) % for MoCap or norm data structure
% mcarrayPlotTimeseries(d, marker, 'dim', dim) % specifying dimensions to be plotted
% mcarrayPlotTimeseries(d, marker, 'timetype', timetype) % axis unit
% mcarrayPlotTimeseries(d, marker, 'plotopt', plotopt) % combined or separate plots
% mcarrayPlotTimeseries(d, marker, 'label', label) % y-axis label 
% mcarrayPlotTimeseries(d, marker, 'names', names) % marker names 
% mcarrayPlotTimeseries(s, segm, 'var', var) % for segm data structure
%
% input parameters
% d/s: MoCap data struct array, norm data struct array, or segm data struct array
% marker: vector containing marker numbers or cell array containing marker names (for MoCap or norm data structure)
% segm: body segment numbers or cell array containing segment names (for segm data structure)
% dim: dimensions to be plotted (for MoCap data structure - default: 1)
% var: variable to be plotted for segment segm (for segm data structure - default: 1)
% timetype: time type used in the plot ('sec' (seconds - default) or 'frame')
% plotopt: plotting option (for MoCap or norm data structure); 'sep' (default) or 'comb':
%   sep: all time series  are plotted in separate subplots
%   comb: all time series will be plotted into the same plot using different colors)
% label: y-axis label (default: no y-axis label). X-axis label is always set, according to timetype
%   (however, for plotting neither x-axis nor y-axis labels: 'label', 0)
% names: if marker names (instead of numbers) are plotted in title and legend (0: numbers (default), 1: names)
% plotMean: plot (overlay) the mean value across the struct array (0 or 1)
% plotStd: plot (overlay) the standard deviation across the struct array (0 or 1)
% showLines: show the individual lines – set to 0 if you want to only display the mean/std curves (0 or 1)
%
%
% output
% Figure.
%
% examples
% mcarrayPlotTimeseries(d, 2) % MoCap or norm data structure, marker 2, dim 1
% mcarrayPlotTimeseries(d, {'Head_FL','Finger_L'}) %marker names instead of numbers (works for segments as well)
% mcarrayPlotTimeseries(d, 1:3, 'dim', 1:3) % markers 1 to 3, dimensions 1 to 3
% mcarrayPlotTimeseries(d, 1:3, 'dim', 3, 'timetype', 'frame') % frames as x axis unit
% mcarrayPlotTimeseries(d, 5, 'dim', 1:3, 'plotopt', 'comb') % all in one plot, different colors per dim
% mcarrayPlotTimeseries(d, 5, 'dim', 1:3, 'plotopt', 'comb', 'label', 'mm') % y-axis label: mm
% mcarrayPlotTimeseries(d, 5, 'dim', 1:3, 'timetype', 'frame', 'label', 0) % no x- axis (and no y-axis) label
% mcarrayPlotTimeseries(d, 5, 'names', 1) % marker names (instead of numbers) plotted in title and legend
% mcarrayPlotTimeseries(s, [3 6 20], 'var', 'angle') % for segm data structure
% mcarrayPlotTimeseries(s, 5:10, 'var', 'eucl', 'timetype', 'frame') % frames as x axis unit
% mcarrayPlotTimeseries(s, [12 14], 'var', 'quat', 'dim', 2, 'plotopt', 'comb') % all in one plot, component 2
% 
%
% 



d=varargin{1};
marker=varargin{2};

dim=[];
var=[];
timetype=[];
plotopt=[];
label=[];
names=[];
plotMean=0;
plotStd=0;
showLines=1;

for k=3:2:length(varargin)
    if strcmp(varargin{k}, 'dim')
        dim=varargin{k+1};
    elseif strcmp(varargin{k}, 'var')
        var=varargin{k+1};
    elseif strcmp(varargin{k}, 'timetype')
        timetype=varargin{k+1};
    elseif strcmp(varargin{k}, 'plotopt')
        plotopt=varargin{k+1};
    elseif strcmp(varargin{k}, 'label')
        label=varargin{k+1};
    elseif strcmp(varargin{k}, 'names')
        names=varargin{k+1};
    elseif strcmpi(varargin{k}, 'plotMean')
        plotMean=varargin{k+1};
    elseif strcmpi(varargin{k}, 'plotStd')
        plotStd=varargin{k+1};
    elseif strcmpi(varargin{k}, 'showLines')
        showLines=varargin{k+1};
    else
        
        str=sprintf('Input argument %s unknown.', varargin{k});
        disp([10, str, 10])
%         [y,fs] = audioread('mcsound.wav');
%         sound(y,fs);
%         return
    end
end

if plotStd
    plotMean = 1;
end

if ~plotMean 
    showLines = 1;
end

%set default values and check for incorrect spelling

if isempty(dim)
    dim=1;
end

if strcmp(d(1).type, 'MoCap data') || strcmp(d(1).type, 'norm data')
    if max(dim)>size(d(1).data,2)/d(1).nMarkers
        disp([10, 'Dimension (dim) value exceeds existing dimensions (no plot created).', 10])
        [y,fs] = audioread('mcsound.wav');
        sound(y,fs);
        return
    end
end

if strcmp(d(1).type, 'segm data')
    if isempty(var)
        disp([10, 'Please specify segment variable (var) to be plotted (no plot created).', 10])
        [y,fs] = audioread('mcsound.wav');
        sound(y,fs);
        return
    end
    if strcmp(var,'angle')
        if max(dim)>1
            disp([10, 'Dimension (dim) value exceeds existing dimensions (no plot created).', 10])
            [y,fs] = audioread('mcsound.wav');
            sound(y,fs);
            return
        end
    elseif strcmp(var,'eucl')
        if max(dim)>3
            disp([10, 'Dimension (dim) value exceeds existing dimensions (no plot created).', 10])
            [y,fs] = audioread('mcsound.wav');
            sound(y,fs);
            return
        end
    elseif strcmp(var,'quat')
        if max(dim)>4
            disp([10, 'Component (dim) value exceeds existing components (no plot created).', 10])
            [y,fs] = audioread('mcsound.wav');
            sound(y,fs);
            return
        end
    end
end

if isempty(timetype)
    timetype='sec';
end
if strcmp(timetype,'sec') || strcmp(timetype,'frame')
else
    timetype='sec';
    disp([10, 'Incorrect spelling of timetype input. Value set to "sec".', 10])
end

if isempty(plotopt)
    plotopt='sep';
end
if strcmp(plotopt,'sep') || strcmp(plotopt,'comb')
else
    plotopt='sep';
    disp([10, 'Incorrect spelling of plotopt input. Value set to "sep".', 10])
end

if isempty(names)
    names=0;
end
if names>1
    disp([10, 'Names input bigger than 1 or in incorrect format, value is set to 1.', 10])
    names=1;
end

if strcmp(d(1).type, 'MoCap data') || strcmp(d(1).type, 'norm data') || strcmp(d(1).type, 'segm data')
else
    disp([10, 'The first input argument should be a variable with MoCap, norm, or segm data structure (no plot created).', 10]);
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    return
end

if nargin<2
    disp([10, 'Please specify data and/or markers (or segments) to be plotted (no plot created).', 10])
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    return
end

if iscell(marker)
    marker=name2number(d(1),marker);
    if isempty(marker) %no plot created if all markers that are given are misspelled
        disp([10, 'Please specify data and/or markers (or segments) to be plotted (no plot created).', 10])
        [y,fs] = audioread('mcsound.wav');
        sound(y,fs);
        return
    end
end

if ischar(marker) %if string is entered
    disp([10, 'Marker number has to be either numbers or cell array with marker names (no plot created).', 10]);
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    return
end
if max(marker)>d(1).nMarkers
    disp([10, 'Marker number out of range (no plot created).', 10]);
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    return
end


p1=marker;
p2=dim;

set(0,'DefaultTextInterpreter','none') %for underscores (and such) in marker names

colors={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

if plotMean %lighter colors when showing mean curves
    colors2 = colors;
    colors  = cellfun(@(x) x*.3+.7,colors,'UniformOutput',false);
end

figure

if isfield(d(1),'type')
    
    for di = 1:length(d)
    
    t = (1:d(di).nFrames)';%-1 taken away as it caused time series to start at 0... [BB 20110301]
    if strcmp(timetype,'sec') 
        t = (t-1)/d(di).freq; 
    end
    if strcmp(d(di).type, 'MoCap data')

        al=1;%amount of lines - for 'comb' plotting
        for k=1:length(p1)
            for m=1:length(p2)
                if strcmp(plotopt, 'sep')
                    subplot(length(p1), length(p2), length(p2)*(k-1)+m),hold on

                    if showLines
                        pl=plot(t, d(di).data(:,3*p1(k)-3+p2(m)));
                        set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    end                   
                    
                    axis([min(t) max(t) -Inf Inf])
                    if names==0
                        title(['Marker ' num2str(p1(k)) ', dim. ' num2str(p2(m))])
                    elseif names==1
                        title(['Marker ' char(d(di).markerName{p1(k)}) ', dim. ' num2str(p2(m))])
                    end
                    
                else
                    hold on
                    if showLines
                        pl=plot(t, d(di).data(:,3*p1(k)-3+p2(m)));
                    
                        set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines

                        if names==0
                            st{al} = ['M. ' num2str(p1(k)) ', dim. ' num2str(p2(m))];
                        elseif names==1
                            st{al} = ['M. ' char(d(di).markerName{p1(k)}) ', dim. ' num2str(p2(m))];
                        end 
                        al=al+1;
                    end
                    
                    if names==0
                        title(['Marker [' num2str(p1) '], dim. [' num2str(p2) ']'])
                    elseif names==1
                        title(['Marker [' num2str(p1) '], dim. [' num2str(p2) ']'])
                    end 
                    axis([min(t) max(t) -Inf Inf])
                    
                    hold on
                end
                if ~isscalar(label) %if label is [] or string, then plot labels
                    if strcmp('sec',timetype)
                        xlabel('seconds')
                    else xlabel('frames')
                    end
                    if ~isempty(label)
                        ylabel(label)
                    end
                end
            end
        end
        if showLines && strcmp(plotopt, 'comb') %plot legend
            leg=legend(st, 'Location', 'EastOutside'); 
        end 
    elseif strcmp(d(di).type, 'norm data')
        
        al=1;%amount of lines - for 'comb' plotting
        %plot(t, d(di).data(:,p1));
        for k=1:length(p1)
            for m=1%:length(p2)
                if strcmp(plotopt, 'sep') 
                    subplot(length(p1), 1, k),hold on
                    if showLines
                        pl = plot(t, d(di).data(:,p1(k)));
                        set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines
                    end
                    
                    axis([min(t) max(t) -Inf Inf])
                    if names==0
                        title(['Marker ' num2str(p1(k)) ', norm data'])
                    elseif names==1
                        title(['Marker ' char(d(di).markerName{p1(k)}) ', norm data'])
                    end
                    
                else
                    hold on
                    if showLines
                        pl=plot(t, d(di).data(:,p1(k))); %FIXBB110102: 'comb' also for norm data
                    
                        axis([min(t) max(t) -Inf Inf])
                        set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines

    %                     title(['Marker [' num2str(p1) '], norm data'])
                        if names==0
                            st{al} = ['M. ' num2str(p1(k))];
                        elseif names==1
                            st{al} = ['M. ' char(d(di).markerName{p1(k)})];
                        end

                        al=al+1;
                    end
                    
                    if names==0
                        title(['Marker [' num2str(p1) '], norm data'])
                    elseif names==1
                        title(['Marker [' num2str(p1) '], norm data'])
                    end
                    hold on
                end
                if ~isscalar(label) %if label is [] or string, then plot labels
                    if strcmp('sec',timetype)
                        xlabel('seconds')
                    else
                        xlabel('frames')
                    end
                    if ~isempty(label)
                        ylabel(label)
                    end
                end
            end
        end
        if showLines && strcmp(plotopt, 'comb') %plot legend
            legend(st, 'Location', 'EastOutside'); 
        end 
    
    elseif strcmp(d(di).type, 'segm data')
        
        tmp=[];
        al=1;%amount of lines - for 'comb' plotting
        for k=1:length(p1)
            tmp = getfield(d(di).segm(p1(k)),var);
            if ~isempty(tmp)
                if strcmp(var, 'angle')
                    if strcmp(plotopt, 'sep')
                        % k
                        subplot(length(p1),1,k),hold on
                        pl = plot(t, tmp);
                        set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines
                        axis([min(t) max(t) -Inf Inf])
                        if names==0
                            title(['Segm. ' num2str(p1(k)) ' - angle'])
                        elseif names==1
                            title(['Segm. ' char(d(di).segmentName{p1(k)}) ' - angle'])
                        end                
                    else
                        hold on
                        pl=plot(t, tmp); %FIXBB201202010: 'comb' also for segm data
                        axis([min(t) max(t) -Inf Inf])
                        set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines
                        if names==0
                            title(['Segm. [' num2str(p1) '] - angle'])
                            st{al} = ['Segm. ' num2str(p1(k))];
                        elseif names==1
                            title(['Segm. [' num2str(p1) '] - angle'])
                            st{al} = ['Segm. ' char(d(di).segmentName{p1(k)})];
                        end 
                        al=al+1;
                        hold on
                    end
                    if ~isscalar(label) %if label is [] or string, then plot labels
                        if strcmp('sec',timetype)
                            xlabel('seconds')
                        else xlabel('frames')
                        end
                        if ~isempty(label)
                            ylabel(label)
                        end
                    end
                elseif strcmp(var, 'eucl')
                    for m=1:length(p2)
                        if strcmp(plotopt, 'sep')
                            subplot(length(p1), length(p2), length(p2)*(k-1)+m),hold on
                            pl = plot(t, tmp(:,p2(m)));
                            set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines
                            if names==0
                                title(['Segm. ' num2str(p1(k)) ', dim. ' num2str(p2(m)) ' - eucl']);
                            elseif names==1
                                title(['Segm. ' char(d(di).segmentName{p1(k)}) ', dim. ' num2str(p2(m)) ' - eucl'])
                            end
                            axis([min(t) max(t) -Inf Inf])
                        else
                            hold on
                            pl=plot(t, tmp(:,p2(m))); %FIXBB201202010: 'comb' also for segm data
                            axis([min(t) max(t) -Inf Inf])
                            set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines
                            if names==0
                                title(['Segm. [' num2str(p1) '], dim. [' num2str(p2) '] - eucl'])
                                st{al} = ['Segm. ' num2str(p1(k)) ', dim. ' num2str(p2(m)) ];
                            elseif names==1
                                title(['Segm. [' num2str(p1) '], dim. [' num2str(p2) '] - eucl'])
                                st{al} = ['Segm. ' char(d(di).segmentName{p1(k)}) ', dim. ' num2str(p2(m)) ];
                            end
                            al=al+1;
                            hold on
                        end
                        if ~isscalar(label) %if label is [] or string, then plot labels
                            if strcmp('sec',timetype)
                                xlabel('seconds')
                            else xlabel('frames')
                            end
                            if ~isempty(label)
                                ylabel(label)
                            end
                        end
                    end
                elseif strcmp(var, 'quat')
                    for m=1:length(p2)
                        if strcmp(plotopt, 'sep')
                            subplot(length(p1), length(p2), length(p2)*(k-1)+m),hold on
                            plot(t, tmp(:,p2(m)))
                            if names==0
                                title(['Segm. ' num2str(p1(k)) ', comp. ' num2str(p2(m)) ' - quat']);
                            elseif names==1
                                title(['Segm. ' char(d(di).segmentName{p1(k)}) ', comp. ' num2str(p2(m)) ' - quat'])
                            end
                            axis([min(t) max(t) -Inf Inf])
                        else
                            hold on
                            pl=plot(t, tmp(:,p2(m))); %FIXBB201202010: 'comb' also for segm data
                            axis([min(t) max(t) -Inf Inf])
                            set(pl,'color',colors{mod(al-1,length(colors))+1}) %al has 4 colors, should start over with blue after 7 lines
                            if names==0
                                title(['Segm. [' num2str(p1) '], comp. [' num2str(p2) '] - quat'])
                                st{al} = ['Segm. ' num2str(p1(k)) ', comp. ' num2str(p2(m)) ];
                            elseif names==1
                                title(['Segm. [' num2str(p1) '], comp. [' num2str(p2) '] - quat'])
                                st{al} = ['Segm. ' char(d(di).segmentName{p1(k)}) ', comp. ' num2str(p2(m)) ];
                            end
                            al=al+1;
                            hold on
                        end
                        if ~isscalar(label) %if label is [] or string, then plot labels
                            if strcmp('sec',timetype)
                                xlabel('seconds')
                            else xlabel('frames')
                            end
                            if ~isempty(label)
                                ylabel(label)
                            end
                        end
                    end
                end
            else
                disp([10, 'No data to be plotted.', 10])
            end
        end
    end
    if showLines && strcmp(plotopt, 'comb') %plot legend
        legend(st, 'Location', 'EastOutside'); 
    end
    end
    
    if plotMean
    al=1;       
    t = (1:min([d.nFrames]))';%-1 taken away as it caused time series to start at 0... [BB 20110301]
    if strcmp(timetype,'sec') 
        t = (t-1)/d(1).freq; 
    end
    if strcmp(d(1).type, 'MoCap data')
        
        dmean = mcarrayMean(d);

        for k=1:length(p1)
            for m=1:length(p2)
                if strcmp(plotopt, 'sep')
                    subplot(length(p1), length(p2), length(p2)*(k-1)+m)

                    pl=plot(t, dmean.data(:,3*p1(k)-3+p2(m)),'k','LineWidth',2);
                    
                else
                    pl=plot(t, dmean.data(:,3*p1(k)-3+p2(m)),'k','LineWidth',2);
                    set(pl,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    if names==0
                        st{al} = ['M. ' num2str(p1(k)) ', dim. ' num2str(p2(m)), ' MEAN'];
                    elseif names==1
                        st{al} = ['M. ' char(d(di).markerName{p1(k)}) ', dim. ' num2str(p2(m)), ' MEAN'];
                    end 
                    al=al+1;
                    hold on
                end
            end
        end
    elseif strcmp(d(1).type, 'norm data')
        
        dmean = mcarrayMean(d);
        
        al=1;%amount of lines - for 'comb' plotting
        %plot(t, dmean.data(:,p1),'k','LineWidth',2);
        for k=1:length(p1)
            for m=1%:length(p2)
                if strcmp(plotopt, 'sep') 
                    subplot(length(p1), 1, k)
                    plot(t, dmean.data(:,p1(k)),'k','LineWidth',2);
                else
                    pl=plot(t, dmean.data(:,p1(k)),'k','LineWidth',2); %FIXBB110102: 'comb' also for norm data
                    set(pl,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    if names==0
                        st{al} = ['M. ' num2str(p1(k)), ' MEAN'];
                    elseif names==1
                        st{al} = ['M. ' char(d(di).markerName{p1(k)}), ' MEAN'];
                    end

                    al=al+1;
                    hold on
                end
            end
        end
    
    elseif strcmp(d(1).type, 'segm data')
        
        disp('plotting of array mean not implemented for segment data')
    end
    if strcmp(plotopt, 'comb') %plot legend
        leg=legend(st, 'Location', 'EastOutside'); 
    end 
                
    end %of if plotmean
    
    
    
    if plotStd
    al=1;    
    if strcmp(d(1).type, 'MoCap data')
        
        dstd = mcarrayStd(d);

        for k=1:length(p1)
            for m=1:length(p2)
                if strcmp(plotopt, 'sep')
                    subplot(length(p1), length(p2), length(p2)*(k-1)+m)

                    pl1 = plot(t, dmean.data(:,3*p1(k)-3+p2(m))-dstd.data(:,3*p1(k)-3+p2(m)),'k','HandleVisibility','off');
                    pl2 = plot(t, dmean.data(:,3*p1(k)-3+p2(m))+dstd.data(:,3*p1(k)-3+p2(m)),'k','HandleVisibility','off');
                    %set(pl1,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    %set(pl2,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                else
                    pl1 = plot(t, dmean.data(:,3*p1(k)-3+p2(m))-dstd.data(:,3*p1(k)-3+p2(m)),'HandleVisibility','off');
                    pl2 = plot(t, dmean.data(:,3*p1(k)-3+p2(m))+dstd.data(:,3*p1(k)-3+p2(m)),'HandleVisibility','off');
                    set(pl1,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    set(pl2,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    if names==0
                        st{al} = ['M. ' num2str(p1(k)) ', dim. ' num2str(p2(m)), ' std'];
                    elseif names==1
                        st{al} = ['M. ' char(d(di).markerName{p1(k)}) ', dim. ' num2str(p2(m)), ' std'];
                    end 
                    al=al+1;
                    hold on
                end
            end
        end
    elseif strcmp(d(1).type, 'norm data')
        
        dstd = mcarrayStd(d);
        
        al=1;%amount of lines - for 'comb' plotting
        %plot(t, dmean.data(:,p1),'k','LineWidth',2);
        for k=1:length(p1)
            for m=1%:length(p2)
                if strcmp(plotopt, 'sep') 
                    subplot(length(p1), 1, k)
                    plot(t, dmean.data(:,p1(k))-dstd.data(:,p1(k)),'HandleVisibility','off');
                    plot(t, dmean.data(:,p1(k))+dstd.data(:,p1(k)),'HandleVisibility','off');
                else
                    pl1 = plot(t, dmean.data(:,p1(k))-dstd.data(:,p1(k)),'HandleVisibility','off');
                    pl2 = plot(t, dmean.data(:,p1(k))+dstd.data(:,p1(k)),'HandleVisibility','off');
                    set(pl1,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    set(pl2,'color',colors2{mod(al-1,length(colors))+1}) %al has 7 colors, should start over with blue after 7 lines
                    
                    if names==0
                        st{al} = ['M. ' num2str(p1(k)), ' std'];
                    elseif names==1
                        st{al} = ['M. ' char(d(di).markerName{p1(k)}), ' std'];
                    end

                    al=al+1;
                    hold on
                end
            end
        end
    
    elseif strcmp(d(1).type, 'segm data')
        
        disp('plotting of array SD not implemented for segment data')
    end
                
    end %of if plotstd
    
    
    
    
else % direct reference to data
    if ~exist('p1')
        p1=1;
    end
    d=d(:);
    plot(d(:,p1));
end
hold off
end


function m=name2number(d, marker)
m=[];


for i=1:length(marker)
    if isfield(d,'markerName') %mocap or norm structure
        for j=1:length(d.markerName)
            if strcmp(char(marker{i}),char(d.markerName{j}))
                m=[m,j];
            end
        end
    elseif isfield(d,'segmentName') %segment structure
        for j=1:length(d.segmentName)
            if strcmp(char(marker{i}),char(d.segmentName{j}))
                m=[m,j];
            end
        end
    end
end

if length(m)<length(marker)
    x=length(marker)-length(m);
    str=sprintf([10 '%d marker(s) not matching with marker names in the given MoCap Structure.', x, 10]);
    disp([10, str, 10])
end

end
