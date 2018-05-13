function warpedData = knMultiDTW(varargin)
% Dynamic time warping function for multiple input signals. All signals are
% warped to the timeline of one of the other input signals.
%
% Use a cell array as input, where each cell contains a time series. 
% Last argument specifies sync array number (Default = 1st)

% Checking if an argument is given to specify which signal to sync to
if isnumeric(varargin{end}) && length(varargin{end}) == 1
    warpIndex = varargin{end}; %index of the time series to which other time series are synced
    varargin(end) = [];
else
    warpIndex = 1; %default: first input signal
end

%converting to cell array if a normal array is used as input signal (rather than cell array)
if iscell(varargin{1})
    varargin = varargin{1};
end

% initializing output matrix
warpedData = zeros(length(varargin),length(varargin{warpIndex}));

%loop through all input signals, warp each one to the reference signal
for i = 1:length(varargin)

    if i == warpIndex
        warpedData(i,:) = varargin{i};
    else
        warpedData(i,:) = knwarp(varargin{warpIndex},varargin{i});
    end

end

end


function b2 = knwarp(a,b)

    %simple pairwise time warping
    [~,ixAB,iyAB] = dtw(a,b);
    
    %initialising output
    b2 = nan(1,length(a));
    
    da = ~logical(diff(iyAB)); %differentiating to check if the index stays the same or if it moves
    db = ~logical(diff(ixAB)); %differentiating to check if the index stays the same or if it moves
    da(end)=0;
    db(end)=0;
   
    
    b2(1)=b(1); %first index of output is equal to first index of input
    i = 2;
    ioffset = 0;   %variable for tracking offsets 

    while i < length(da)+1
        
        if da(i) %if the time index of signal a stays the same
            if iyAB(i)-ioffset > 1
                b2(iyAB(i)-ioffset)= b2(iyAB(i)-ioffset-1); %set value equal to previous index
                ioffset=ioffset-1;                          %decrease offset value
            end
            i = i+1;
            
        elseif db(i) %if the time index of signal b stays the same
            xlow = i;
            ii = 0; %find out how long it stays the same
            while db(i + ii)
                ii = ii+1;
            end
            xhigh = ii+i;
            
            b2(iyAB(i)-ioffset) = mean(b(iyAB(xlow:xhigh))); %set output to the mean value within this period
            ioffset = ii+ioffset;                            %increase offset value
            i = i+ii+1;
        else        %the time indecies of both signals increases
            b2(iyAB(i)-ioffset) = b(iyAB(i));
            i = i+1;
        end
    
        
    end
    

end
