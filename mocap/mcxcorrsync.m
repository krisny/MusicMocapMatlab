function dout = mcxcorrsync(d,syncper)
%
% Syncronize streams by cross correlation
%
% d: mocap structure
% syncper: start and end time (in seconds) of sync movement
%
% dout: new mocap struct where all markers are aligned to the first marker
%
%


%todo: option for calculating timeder before sync
%todo: optional reference marker (currently only using first marker)


dout = d;
dout.data = NaN(size(dout.data,1),size(dout.data,2));

syncper = d.freq*syncper; %converting syncper to frame numbers



if ~strcmp(d.type,'norm data')
    %todo: option for not calculating norm before sync
    d2 = mcnorm(d);
    dout.data(:,1) = d.data(:,1);
else
    d2 = d;
    dout.data(:,1:3) = d.data(:,1:3);
end




%creating NaN matrix (ensures that empty frames becomes NaN after aligning)



for i = 2:d2.nMarkers
    
    %find lag
    [r,lags] = xcorr(d2.data(syncper(1):syncper(2),1),d2.data(syncper(1):syncper(2),i));
    [~,l ] = max(r);
    
    %%%keeping the plotting code for now for debugging%%%
    %subplot(4,5,i)
    %plot(lags,r)
    %hold all
    %scatter(lags(l),r(l))
    %hold off
    %title(['ab lag: ' num2str(lags(l)) ' samples. corr: ' num2str(r(l)) ''])
    
    
    %timeshift by estimated lag
    if strcmp(d.type,'norm data')
        dout.data(trimtorange(d2.nFrames,lags(l)),i) = d.data(trimtorange(d2.nFrames,-lags(l)),i);
    else
        xyz = i*3-[2 1 0];
        dout.data(trimtorange(d2.nFrames,lags(l)),xyz) = d.data(trimtorange(d2.nFrames,-lags(l)),xyz);
    end
end



end


function datarange = trimtorange(nFrames,lag)

datarange = (1:nFrames)+lag;

datarange(datarange < 1) = [];
datarange(datarange > nFrames) = [];

end

