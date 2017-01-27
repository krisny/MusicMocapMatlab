
function mcrtanimate(d,p,speedlim)
%%
% Function for real time 3D animation of mocap data, with sync to audio
% 
% mcrtanimate(d)
% mcrtanimate(d,p)
% mcrtanimate(d,p,speedlim)
%
% 
% d is a MoCap struct
%
% p is a animpar struct, where the .conn field is used for making bones in
% the animation
%
% speedlim is a time limitation in seconds for the streaming of sensor data 
% to Max (default = 0.05)
%
%
% Use interface in Max to control the animation:
% Play, stop, tempo, loop on/off, loop marks, 
%
% If audio sync is desired, one extra field must be added to the mocap
% struct:
% 
% d.other.soundFile = 'filename' (without path information)
% By default, the path is chosen to be a folder named 'sounds' in the
% current working folder. If another folder is desired, this should be
% chosen in the Max GUI.
%
% For now, the mcrtanimate.maxpat must be opened manually in order for the
% function to work. Will try to make this automated later.
%
%
% by 
% Kristian Nymoen
% University of Oslo 
% fourMs - Music, Mind, Motion, Machines
%
% Addon to the Motion Capture Toolbox, University of Jyvaskyla, Finland

if nargin < 2
    p.conn = [];
    speedlim = 0.05; %speed limit in seconds
elseif nargin < 3
    speedlim = 0.05; %speed limit in seconds
end

global playon
global currentFrame
global startPoint
global endPoint
global loop
global timestretch
global closing
 
playon = [0 0];
closing = 0;
%speedlim = 0.1; %speed limit in seconds

UDPin = udp('127.0.0.1',7703);
UDPout = udp('127.0.0.1',7701);


UDPin.LocalPort = 7702;
UDPout.OutputBufferSize = 2048;         %temporary workaround to be able to stream more data...
UDPout.OutputDatagramPacketSize = 2048; %temporary workaround to be able to stream more data...

%% 
% Setup a figure window and define a callback function for close operation
figureHandle = figure('NumberTitle','off',...
    'Name','mcrtanimate - close this fig to close UDP connection',...
    'Color',[0 0 0],'UserData',0,...
    'CloseRequestFcn',{@localCloseFigure,UDPin,UDPout});

%%
% Define a callback function to be executed when a datagram is received
UDPin.DatagramReceivedFcn = {@readFromMaxGUI,figureHandle};

%% 
% Open the interface object
fopen(UDPout);
fopen(UDPin);

pause(1);

%%
% Send initial data information
fwrite(UDPout,sprintf('/info/fileName %s',d.filename));
if isfield(d,'other')
    if isfield(d.other,'soundFile')
        fwrite(UDPout,sprintf('/info/soundFile %s',d.other.soundFile));
    end
end
fwrite(UDPout,sprintf('/info/nMarkers %d',d.nMarkers));
pause(1);
fwrite(UDPout,sprintf('/info/freq %d',d.freq));
fwrite(UDPout,sprintf('/info/length %f',d.nFrames/d.freq));
fwrite(UDPout,sprintf('/info/currentPath %s%s',pwd,'/sounds'));
fwrite(UDPout,sprintf('/info/bones %s',num2str(reshape(p.conn',1,numel(p.conn)))));

fwrite(UDPout,sprintf('/info/nFrames %d',d.nFrames));


%%
%give it som time..
pause(1);




%%
%infinite loop as long as closing == 0
while ~closing 

    %Start from startpoint
    currentFrame = startPoint;

    
    while currentFrame <= endPoint 
        if ~playon(1);
            playon(2) = 0;
            waitfor(figureHandle,'UserData')
        end
        tic

        if ~closing
            %output currentFrame
            fwrite(UDPout,sprintf('/info/currentFrame %d',currentFrame));

            %output mocapdata
            fwrite(UDPout,sprintf('/animation/data %s',num2str(d.data(currentFrame,:))));


            pause(timestretch*1/d.freq-toc*1.2) %necessary to scale toc by 1.2 for some strange reason...
            
            
            if toc > speedlim
                currentFrame = currentFrame+1;
            else
                %do something to skip frames...
                currentFrame = currentFrame+round((speedlim-toc)*d.freq/timestretch);%(speedlim-toc);
                
                pause(speedlim-toc)
            end

        end
        
    end
 
if loop == 0 && ~closing
    %if looping is not activated, then stop
    playon = [0 0];
    fwrite(UDPout,sprintf('/info/play %d',0));
end

end



%% callback when datagram is received from Max
function readFromMaxGUI(UDPin,~,figureHandle)

%%
%need clobal variables to be able to set them in the callback function
global playon
global currentFrame
global startPoint
global endPoint
global loop
global timestretch
global closing

%disp('callback')

%% 
% read the variables
% string: play, currentFrame, startPoint, endPoint, loop, tempo, closing
data = str2num(fscanf(UDPin));

%% 
% set the varables
playon(1) = data(1);
currentFrame = data(2);
startPoint = data(3);
endPoint = data(4);
loop = data(5);
timestretch = data(6);
closing = data(7);
set(figureHandle,'UserData',currentFrame);

if playon(1) && ~playon(2)
    set(figureHandle,'UserData',currentFrame+1);
    playon(2) = 1;
    if currentFrame == endPoint
        currentFrame = startPoint;
    end
%    uiresume
end

if closing
    currentFrame = endPoint+1;
    set(figureHandle,'UserData',currentFrame+1);
    close(figureHandle)
    return
end




%% Implement the close figure callback - when the figure closes, the UDP connections are closed
function localCloseFigure(figureHandle,~,UDPin,UDPout)

global currentFrame
global endPoint
global closing

fwrite(UDPout,sprintf('/info/matlabClose bang'));

closing = 1;
currentFrame = endPoint+1;

fclose(UDPin);
fclose(UDPout);
delete(UDPin);
fclose(UDPout);
clear UDPin UDPout;

delete(figureHandle);
