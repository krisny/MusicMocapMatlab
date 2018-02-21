function k = pitchExtract(varargin)
%Extracting pitch envelope from a monophone audio source.
%
% Usage:
% pitchExtract(d)
% d may be one of the following:
%  - a string indicating a file name
%  - a miraudio object (MIR toolbox)
%  - a variable containing audio samples at 44100 Hz
%  - the string 'record' records and analyses 10 seconds of audio from the microphone
%
% pitchExtract(d,varargin)
%
% Various arguments can be given to adjust the function.
% pitchExtract(d,'playsound',1) - plays the sound file along with analysis result (Default)
% pitchExtract(d,'playsound',0) - does not play anything
%
% pitchExtract(d,'limits',[t1 t2]) - adjusts the time limits for playback of audio and plotting.
%	t1 and t2 in seconds. Plays and plots entire sound file by default.
%
% pitchExtract(d,'fundamental',f) - adjusts the y-axis in the plotting function,
%	showing just/pure intonation intervals rather than the actual frequency.
%
% pitchExtract(d,'volume',f) - adjusts the volume of the analysis sine tone
%   during playback. Default value is 0.25.
%
% Other arguments that can be changed are various parameters and thresholds:
%	'windowsize'        - size of window used in autocorrelation function
%                         (default 1/10 of sampling rate)
%	'hop'               - hop size used in autocorrelation function
%                         (default 0.1)
%	'freqlim'           - frequency limits given to the pitch extraction
%                         function (default [180 750])
%	'maxOctaveJump'     - Threshold for correcting octave jumps in pitch
%                         extraction algorithm (default 20 frames).
%	'dThresh1'          - Threshold for segmentation algorithm (default 1)
%	'dThresh2'          - Threshold for segmentation algorithm (default 5)
%	'filter1'           - Filter length for segmentation algorithm
%                         (default 15)
%	'filter2'           - Filter length for segmentation algorithm
%                         (default 21)
%	'toneLengthThresh'  - Threshold for segmentation algorithm (default 10)
%	'rmsThresh'         - Thresholding RMS envelope in segmentation
%                         algorithm (default 0.005)
%	'corrThresh'        - Thresholding autocorrelation coefficient for
%                         segmentation algorithm (default 0.7)
%
%
% By Kristian Nymoen, 2017


d=varargin{1};

if ischar(d)
    if strcmp(d,'record')
        recObj = audiorecorder(44100,16,1);
        disp('Now recording');
        recordblocking(recObj,10)
        disp('Ended recording. Analysing...');
        lyd = getaudiodata(recObj);
        k.sr = 44100;
    else
        [lyd,k.sr] = audioread(d);
    end
elseif strcmp(mirtype(d),'miraudio')
    k.sr = uncell(get(d,'Sampling'));
    lyd = mirgetdata(d);
elseif isnumeric(d)
    k.sr = 44100;
    lyd = d;
    disp('Treating input as audio signal with sampling rate 44100 Hz')
end

k.audio = lyd;

if size(lyd,2) == 2 %if stereo
    k.audio = mean(lyd,2); %make mono
elseif size(lyd,1) == 2 %if stereo
    k.audio = mean(lyd,1)'; %make mono
end



v.plotlim = [1/k.sr length(lyd)/k.sr];
v.fundamental = 1;
v.medianFilterLength1 = 15;
v.derivTreshold1 = 1;
v.medianFilterLength2 = 21;
v.derivTreshold2 = 5;
v.toneLengthTreshold = 10;
v.rmsThreshold = 0.005;
v.ws = round(k.sr/10);
v.hop = 0.1;
v.minf = 180;
v.maxf = 750;
v.playsound = 1;
v.corrThreshold = 0.7;
v.maxOctaveJump = 20;
v.sinusvolume = 0.25;

for argnr=2:2:length(varargin)
    if strcmpi(varargin{argnr}, 'limits')
        v.plotlim=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'fundamental')
        v.fundamental=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'dThresh1')
        v.derivTreshold1=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'dThresh2')
        v.derivTreshold2=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'filter1')
        v.medianFilterLength1=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'filter2')
        v.medianFilterLength2=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'toneLengthThresh')
        v.toneLengthTreshold=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'rmsThresh')
        v.rmsThreshold=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'windowsize')
        v.ws=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'hop')
        v.hop=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'freqlim')
        v.minf=varargin{argnr+1}(1);
        v.maxf=varargin{argnr+1}(2);
    elseif strcmpi(varargin{argnr}, 'playsound')
        v.playsound=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'corrThresh')
        v.corrThreshold=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'maxOctaveJump')
        v.maxOctaveJump=varargin{argnr+1};
    elseif strcmpi(varargin{argnr}, 'volume')
        v.sinusvolume=varargin{argnr+1};
    else
        str=sprintf('Input argument %s unknown.', varargin{argnr});
        disp([10, str, 10])
    end
end


[k.corrs,k.pitches,k.steps,k.corrcoeff,k.rms]= windowedAutoCorr(k,v);



%%


[k.mtones,k.durations] = segmentAndPlotPitches(k,v);



%%

playAnalysedAudio(k,v)



%%
%savepitches2(k.mtones',k.durations,['pitchstep' sanger(i).name(1:2) 'midi.txt']);


k.mtones = k.mtones';

if v.fundamental > 1
    k.pitchesCents = 1200 * log2(k.pitches/v.fundamental);
    k.mtonesCents = 1200 * log2(k.mtones/v.fundamental);
end

end



function [mtones,durations] = segmentAndPlotPitches(x,v)

pitches = x.pitches;
steps = x.steps;
rmss = x.rms;

if nargin < 2
    v.x = [];
end

fundamental = 1;
medianFilterLength1 = 15;
derivTreshold1 = 1;
medianFilterLength2 = 21;
derivTreshold2 = 5;
toneLengthTreshold = 10;
rmsThreshold = 0.005;

if isfield(v,'fundamental')
    fundamental = v.fundamental;
end
if isfield(v,'medianFilterLength1')
    medianFilterLength1 = v.medianFilterLength1;
end
if isfield(v,'derivTreshold1')
    derivTreshold1 = v.derivTreshold1;
end
if isfield(v,'medianFilterLength2')
    medianFilterLength2 = v.medianFilterLength2;
end
if isfield(v,'derivTreshold2')
    derivTreshold2 = v.derivTreshold2;
end
if isfield(v,'toneLengthTreshold')
    toneLengthTreshold = v.toneLengthTreshold;
end
if isfield(v,'rmsThreshold')
    rmsThreshold = v.rmsThreshold;
end

%if RMS energy is really low, penaltize the next frame. (prevents system 
%from trying to extract pitches too soon after a silent period)
for i = 1:length(rmss)-1
    if rmss(i) < rmsThreshold*0.2
        rmss(i+1) = rmss(i+1)*0.1;
    end
end

%SEGMENTATION CODE:

%Filter the pitch curve before segmentation. Removes vibrato.
smoothpitches = medfilt1(pitches,medianFilterLength1,'omitnan','truncate');
smoothpitches(isnan(pitches)) = nan;

%Derivate the pitch curve. If the pitch derivative is larger than
%threshold, remove pitch information (in effect: segment curve)
smoothderiv = [0 diff(smoothpitches)];
smoothpitches(abs(smoothderiv)>derivTreshold1) = nan;

%If the RMS energy is lower than threshold, remove pitch information
smoothpitches(rmss < rmsThreshold) = nan;

%Another round of filtering, derivation and thresholding.
smoothpitches = medfilt1(smoothpitches,medianFilterLength2,'omitnan','truncate');
smoothpitches(rmss < rmsThreshold) = nan;
smoothderiv = [0 diff(smoothpitches)];
smoothpitches(abs(smoothderiv)>derivTreshold2) = nan;



j = 1;
l = [0 0];

%Identify tones. 
% If the pitch curve is nan, and the following non-nan period is longer
% than toneLengthThreshold, then identify a new tone.
for i = 1:length(smoothpitches)
    if isnan(smoothpitches(i))
        if l(2) > 0 && l(1)-l(2) > toneLengthTreshold
            tones{j} = smoothpitches(l(2)+1:l(1));
            durations(j,1:2) = steps([l(2)+1,l(1)]);
            j = j+1;
        end
        l(2) = i;
    end
    l(1) = i;
end

%pitch of each tone is the median of all underlying pitches
for i = 1:length(tones)
    mtones(i) = median(tones{i});%/fundamental;
end


%PLOT CURVE

clf



if fundamental ~= 1 %if a fundamental is specified, then use cents as unit
    plot(steps,1200 * log2(pitches/fundamental));%/fundamental) %raw pitch curve
    hold all
    plot(steps,1200 * log2(smoothpitches/fundamental));%/fundamental) %pitch segmentation/analysis curve
    ylabel('Cents')
    line(durations',[1200 * log2(mtones/fundamental);1200 * log2(mtones/fundamental)],'Color','k','LineWidth',3) %segments

else
    plot(steps,pitches); %raw pitch curve
    hold all
    plot(steps,smoothpitches); %pitch segmentation/analysis curve
    ylabel('Frequencies')
    line(durations',[mtones;mtones],'Color','k','LineWidth',3) %segments
    set(gca,'yScale','log')
end

xlabel('Time')


grid on

if isfield(v,'plotlim')
    xlim(v.plotlim)
end

    
end



function [corrs,pitches,steps,acc,rmss] = windowedAutoCorr(k,v)


a = k.audio;
durationinseconds = length(k.audio)/k.sr;


ws = v.ws;
hop = v.hop;
minf = v.minf;
maxf = v.maxf;
corrThreshold = v.corrThreshold;
maxOctaveJump = v.maxOctaveJump;


% subtracting mean:
a=a-mean(a);



%compute maxlag (limit for xcorr function) from sr and lower frequency limit
maxlag=round(1/minf*k.sr);


%initialize correlation matrix
windowStarts = 1:floor(ws*hop):(length(a)-ws);
corrs = zeros(maxlag*2+1,length(windowStarts));
%corrs = [];


%compute autocorrelations in all window x
j=1;
for i = windowStarts
    windowsignal = a(i:(i+ws-1)).*tukeywin(ws,0.7);
    [corrs(:,j),~] = xcorr(windowsignal,windowsignal,'coeff',maxlag);    
    
    %tmp solution for calculating rms:
    rmss(j)=rms(smooth(windowsignal,5));
        
    %increase index
    j=j+1;
end


for i = 1:size(corrs,2)
    [pitches(i), acc(i)] = findautocorrpeak(corrs(:,i),1/maxf*k.sr);%mod(i,3)==0);
    [pitchesUp(i), pitchesDown(i)] = findautocorrpeakScaled(corrs(:,i),1/maxf*k.sr);%mod(i,3)==0);
end


pitches = 1./(pitches/k.sr);
pitchesDown = 1./(pitchesDown/k.sr);
pitchesUp = 1./(pitchesUp/k.sr);



%preventing octave jumps ... could be improved

%if a single frame jumps an octave, change it:
for i = 2:length(pitches)-1
    if abs(pitches(i-1)-pitches(i+1)) < 80
        if  abs((pitches(i-1)+pitches(i+1))/2 - pitches(i)*0.5) < 80
            pitches(i) = pitches(i)*0.5;
        elseif abs((pitches(i-1)+pitches(i+1))/2 - pitches(i)*2) < 80
            pitches(i) = pitches(i)*2;
        end
    end
end



%if a section shorter than maxOctaveJump changes octave, then change it:
for i = 2:length(pitches)-1
    if abs(pitches(i-1)-pitches(i)) > min(pitches(i-1),pitches(i))-50
        j = i+1;
        while j > i && j < i+maxOctaveJump && j<length(pitches)
            if abs(pitches(j+1)-pitches(j)) > min(pitches(j+1),pitches(j))-50
                a = pitches(i:j)*2;
                b = pitches(i:j)*0.5;
                if range([a,pitches(j+1),pitches(i-1)]) < range([b,pitches(j+1),pitches(i-1)])
                    pitches(i:j)=a;
                else
                    pitches(i:j)=b;
                end
                j = 0;
            else
                j = j+1;
            end
        end
    end
end


%remove pitch values outside legal range.
pitches(pitches < minf) = nan;
pitches(pitches > maxf) = nan;
pitchesDown(pitchesDown < minf) = nan;
pitchesDown(pitchesDown > maxf) = nan;
pitchesUp(pitchesUp < minf) = nan;
pitchesUp(pitchesUp > maxf) = nan;

%remove pitch values with correlation coefficients  below threshold
pitches(acc < corrThreshold) = nan;
pitchesDown(acc < corrThreshold) = nan;
pitchesUp(acc < corrThreshold) = nan;

%expand all nan-fields by one in each direction
%initializing 3 temporary variables for pitches, up and down. 
tmppitches = pitches;
tpd=pitchesDown;
tpu=pitchesUp;
for i = 2:length(pitches)-1
    if isnan(pitches(i))
        tmppitches(i-1) = NaN;
        tmppitches(i+1) = NaN;
    end
    if isnan(pitchesDown(i))
        tpd(i-1) = NaN;
        tpd(i+1) = NaN;
    end
    if isnan(pitchesUp(i))
        tpu(i-1) = NaN;
        tpu(i+1) = NaN;
    end
end

pitches = tmppitches;

pitchesDown = tpd;
pitchesUp= tpu;

for i =1:length(pitches)
    if isnan(pitches(i))
        if ~isnan(pitchesUp(i))
            pitches(i)=pitchesUp(i);
        elseif ~isnan(pitchesDown(i))
            pitches(i)=pitchesDown(i);
        end
    end
end


steps=0:((durationinseconds)/(size(corrs,2)-1)):durationinseconds;


end



function [acp, acc] = findautocorrpeak(x,minl)

minl = round(minl);

i0 = (length(x)+1)/2;
i = i0+floor(minl/2);
corrNegative = 0;
localMinimum = 0;

%move pointer to place with negative correlation value
if min(x) < 0
    
    while corrNegative == 0 && i < length(x)
        if x(i) < 0
            corrNegative = 1;
        else
            i = i+1;
        end
    end
    
    %find a local minimum
    while localMinimum == 0  && i < length(x)
        asdf = min([i+floor(minl/2),length(x)]);
        if x(i) == min(x(i:asdf))
            localMinimum = 1;
        else
            i = i+1;
        end
    end
    
    if i - i0 < minl
        i = minl+i0;
    end
        
    if length(x)-2 < i
        acc = nan;
        I = nan;
        acp = nan;
    else
        upsX = interp1(i:length(x),x(i:length(x)),i:0.1:length(x),'spline');        
        [acc, I] = max(upsX);
        acp = I/(length(upsX)/(length(x)-i+1)) + i - i0;
    end
else
    acc = nan;
    acp = nan;
end


%end
end



function [acp1, acp2] = findautocorrpeakScaled(x,minl)

minl = round(minl);

i0 = (length(x)+1)/2;
i = i0+floor(minl/2);
corrNegative = 0;
localMinimum = 0;

%move pointer to place with negative correlation value
if min(x) < 0
    
    while corrNegative == 0 && i < length(x)
        if x(i) < 0
            corrNegative = 1;
        else
            i = i+1;
        end
    end
    
    %find a local minimum
    while localMinimum == 0  && i < length(x)
        asdf = min([i+floor(minl/2),length(x)]);
        if x(i) == min(x(i:asdf))
            localMinimum = 1;
        else
            i = i+1;
        end
    end
    
    if i - i0 < minl
        i = minl+i0;
    end
        
    if length(x)-2 < i
        I = nan;
        acp1 = nan;
        acp2 = nan;
    else
        
        
        upsX = interp1(i:length(x),x(i:length(x)),i:0.1:length(x),'spline');
        
        [~, I] = max(upsX./((1:length(upsX)).^0.025));   %penalizing higher index values..
        [~, I2] = max(upsX./((1:length(upsX)).^-0.025)); %penalizing lower index values..
        
        
        acp1 = I/(length(upsX)/(length(x)-i+1)) + i - i0;
        acp2 = I2/(length(upsX)/(length(x)-i+1)) + i - i0;
    end
else
    acp1 = nan;
    acp2 = nan;
end


end



function playAnalysedAudio(k,v)

if v.playsound
    
    d = k.durations([k.durations(:,2) > v.plotlim(1) & k.durations(:,1) < v.plotlim(2)]',1:2);
    m = k.mtones(1,[k.durations(:,2) > v.plotlim(1) & k.durations(:,1) < v.plotlim(2)]')';
    
    sinsignal = zeros(size(k.audio(round(v.plotlim(1)*44100):round(v.plotlim(2)*44100))));
    j = 1;
    
    if d(1,1) < v.plotlim(1)
        sinsignal(1:length(0:1/44100:(d(1,2)-v.plotlim))) = sin((0:1/44100:(d(1,2)-v.plotlim))*2*pi*m(1));
        j=j+1;
    end
    
    while j <= length(m)
        sinsignal(round((d(j,1)-v.plotlim)*44100):round((d(j,1)-v.plotlim)*44100)+length((0:1/44100:(d(j,2)-d(j,1))))-1) = sin((0:1/44100:(d(j,2)-d(j,1)))*2*pi*m(j));
        j = j+1;
    end
    
    sinsignal = smooth(sinsignal,30)*v.sinusvolume;
    
    if length(length(k.audio(round(v.plotlim(1)*44100):round(v.plotlim(2)*44100)))) < length(sinsignal)
        sinsignal(1+length(k.audio(round(v.plotlim(1)*44100):round(v.plotlim(2)*44100))):end) = [];
    end
    
    
    %sound(sinsignal,44100)
    %%
    sound([...
        k.audio(round(v.plotlim(1)*44100):round(v.plotlim(2)*44100)),...
        sinsignal],...
        44100)
    
end

end