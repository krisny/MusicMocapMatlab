function warpedData = knMultiDTW(varargin)
%cell array as input, each cell contains a time series. Add sync array
%number as an optional argument

if isnumeric(varargin{end}) && length(varargin{end}) == 1
    warpIndex = varargin{end}; %index of the time series to which other time series are synced
    varargin(end) = [];
else
    warpIndex = 1; 
end

%converting to cell array in case someone inputs a normal array
if iscell(varargin{1})
    varargin = varargin{1};
end


warpedData = zeros(length(varargin),length(varargin{warpIndex}));

for i = 1:length(varargin)
    
    if i == warpIndex
        warpedData(i,:) = varargin{i};
    else
        warpedData(i,:) = knwarp(varargin{warpIndex},varargin{i});
    end

end

end


function b2 = knwarp(a,b)


    [~,ixAB,iyAB] = dtw(a,b);
    

    b2 = ones(1,length(a))*3;
    
    da = ~logical(diff(iyAB));
    db = ~logical(diff(ixAB));
    da(end)=0;
    db(end)=0;
   
    
    b2(1)=b(1);
    i = 2;
    negi = 0;

    while i < length(da)+1
    
        if da(i)
            if iyAB(i)-negi > 1
            b2(iyAB(i)-negi)= b2(iyAB(i)-negi-1);
            negi=negi-1;
            end
            i = i+1;
            
        elseif db(i)
            xlow = i;
            ii = 0;
            while db(i + ii)
                ii = ii+1;
            end
            xhigh = ii+i;
            
            b2(iyAB(i)-negi) = mean(b(iyAB(xlow:xhigh))); 
            negi = ii+negi;
            i = i+ii+1;
        else
            b2(iyAB(i)-negi) = b(iyAB(i));
            i = i+1;
        end
    
        
    end
    

end
