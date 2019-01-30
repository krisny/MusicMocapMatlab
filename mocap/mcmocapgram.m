function h = mcmocapgram(varargin)
% Plots mocapgram (shows positions of a large number of markers as projection 
% onto a colorspace).
%
% syntax
% h = mcmocapgram(d);
% mcmocapgram(d);
% mcmocapgram(h,d);
% mcmocapgram(d,timetype);
% mcmocapgram(h,d,timetype);
% mcmocapgram(h,d,maxscale,timetype);
%
% input parameters
% d: MoCap data structure, Norm structure
% timetype: time type used in the plot ('sec' (default) or 'frame')
% maxscale: input values will be normalised between 0 and this value. If
%           not given. Input values will be normalised between min and max
%           for each data column (X, Y, Z for each marker individually).
% h: axes handle (preferred) or figure handle
%
% output
% h: figure handle
%
% examples
% mcmocapgram(d,'frame');
% mcmocapgram(h,d,1,'frame');
% mcmocapgram(h,d,1000);
% mcmocapgram(d,1);
% h = mcmocapgram(d);
%
% Script developed by Kristian Nymoen, University of Oslo, Norway
% ? Part of the Motion Capture Toolbox, Copyright ?2008,
% University of Jyvaskyla, Finland

if ishandle(varargin{1})
    h = varargin{1};
    d = varargin{2};
    if length(varargin) > 2
        if isnumeric(varargin{3})
            maxscale = varargin{3};
        end
    end
else
    d = varargin{1};
    if length(varargin) > 1
        if isnumeric(varargin{2})
            maxscale = varargin{2};
        end
    end
end


if strcmp(d.type, 'MoCap data') || strcmp(d.type, 'norm data')
    if ischar(varargin{end})
        timetype=varargin{end};
    else
        timetype = 'sec';
    end
else
    disp([10, 'The first input argument should be a variable with MoCap data structure.', 10]);
    return
end

r = zeros(d.nMarkers,d.nFrames);
g = zeros(d.nMarkers,d.nFrames);
b = zeros(d.nMarkers,d.nFrames);

if strcmp(d.type, 'MoCap data')
    for i = 1:d.nMarkers
        for j = 1:d.nFrames
            r(i,j) = d.data(j,i*3-2);
            g(i,j) = d.data(j,i*3-1);
            b(i,j) = d.data(j,i*3);
        end
    end
end

if strcmp(d.type, 'norm data')
    for i = 1:d.nMarkers
        for j = 1:d.nFrames
            r(i,j) = d.data(j,i);
            g(i,j) = d.data(j,i);
            b(i,j) = d.data(j,i);
        end
    end
end

if exist('maxscale','var')
    mr = maxscale;
    mg = maxscale;
    mb = maxscale;
    minr = 0;
    ming = 0;
    minb = 0;
end


for i = 1:d.nMarkers
    if ~exist('maxscale','var')
        mr = max(r(i,:));
        mg = max(g(i,:));
        mb = max(b(i,:)); 
        r(i,:) = r(i,:)-min(r(i,:));
        g(i,:) = g(i,:)-min(g(i,:));
        b(i,:) = b(i,:)-min(b(i,:));
    end
    r(i,:) = r(i,:)./mr;
    g(i,:) = g(i,:)./mg;
    b(i,:) = b(i,:)./mb;
end
rgb(:,:,1) = r;
rgb(:,:,2) = g;
rgb(:,:,3) = b;


if exist('h','var')
    try 
        set(0, 'currentfigure', h);
        h2 = axes(h);
        if strcmp(timetype,'frame')
            image([0 d.nFrames],[1 d.nMarkers],rgb)
            xlabel('time (frames)')
        else
            image([0 d.nFrames/d.freq],[1 d.nMarkers],rgb)
            xlabel('time (s)')
        end

    catch
        h2=h;
        if strcmp(timetype,'frame')
            image(h2,[0 d.nFrames],[1 d.nMarkers],rgb)
            xlabel('time (frames)')
        else
            image(h2,[0 d.nFrames/d.freq],[1 d.nMarkers],rgb)
            xlabel('time (s)')
        end
        
    end
else
    h = figure;figure(gcf);
    h2 = axes(h);
    if strcmp(timetype,'frame')
        image([0 d.nFrames],[1 d.nMarkers],rgb)
        xlabel('time (frames)')
    else
        image([0 d.nFrames/d.freq],[1 d.nMarkers],rgb)
        xlabel('time (s)')
    end

end


ylabel('marker')

if ~isempty(d.markerName)
    if ischar(d.markerName{1})
        set(h2,'YTick',1:d.nMarkers,'YTickLabel',d.markerName)
    elseif iscell(d.markerName{1})
        set(h2,'YTick',1:d.nMarkers,'YTickLabel',[d.markerName{:}])
    else
        disp('unknown markerName type')
    end
end





if nargout == 0
    clear h
end    



end
