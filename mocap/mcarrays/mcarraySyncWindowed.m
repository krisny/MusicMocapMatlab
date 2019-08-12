function data = mcarraySyncWindowed(data,masterIndex,wl,hop,maxlag,td)

%
% Function that synchronises a dataset of more or less similar data. e.g. multiple examples of the same choreography
%
% mcarraySyncWindowed(data,masterIndex,windowLength,hopSize,maxLag,timeDerivative)
%
% data: An array of mocap structs 
%
% masterIndex: The index within the dataset to which the other recordings
% will be aligned. If set to zero, the dataset is synced to the mean across all recordings. Default 1.
%
% windowLength: the length (in frames) of the moving data window (default 4000)
% hopSize hop size in frames for the moving data window (default 200)
% 
% maxLag: the maximum lag considered in the cross correlation.
%
% timeDerivative: time derivative order. Use if syncing with e.g. velocity
% or acceleration. (Default 0 - sync using position data)
%
% Sync entire dataset to one entry in the array or to the mean of the array
% Sync strategy: 
%   Calculate windowed cross correlation for all markers, 
%   Window size 4000 frames by default, hop 200 frames by default.
%   Pick correlation peak within +/- 400 (default) frames
%   Median across all correlation peaks
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%




    if nargin <= 2
        masterIndex = 1; 
    end

    if nargin <= 3
        wl = 4000; 
    end
    
    if nargin <= 4
        hop = 200; %hop length        
    end
    
    if nargin < 5
        maxlag = 400;
    end

    if nargin < 6
        td = 0;
    end
    
    padframes = maxlag;
    
    %isolate the master recording, to which the other recordigs are aligned
    if masterIndex == 0
        xcm = mcarrayMean(data);
        if td
            xcm = mctimeder(xcm,td);
        end
    else
        if td
            xcm = mctimeder(data(masterIndex),td);
            xcm = mcfillgaps(xcm);
        else
            xcm = mcfillgaps(data(masterIndex));
        end
        data(masterIndex).syncpoint = 0;
    end

    xcm = mcsmoothen(xcm,15);
    xcm.data = [nan(padframes, xcm.nMarkers*3); xcm.data; nan(padframes, xcm.nMarkers*3)];
    xcm.nFrames = xcm.nFrames+padframes*2;
    
    
    %extract the indecies for the rest of the dataset (excluding master)
    toSync = setdiff(1:length(data),masterIndex);

    for i = toSync    
        
        if td
            xci = mctimeder(data(i),td); 
            xci = mcfillgaps(xci);
        else
            xci = mcfillgaps(data(i));
        end
        
        
        xci = mcsmoothen(xci,15);        
        xci.data = [nan(padframes, xci.nMarkers*3); xci.data; nan(padframes, xci.nMarkers*3)];
        xci.nFrames = xci.nFrames+padframes*2;
        
        syncCandidates = [];
        lagConfidences = [];
        
        for k = padframes:hop:(xcm.nFrames-padframes-wl) %loop through overlapping time windows
            
            allrs = nan(xci.nMarkers*3,maxlag*2+1);
            
            
            for j = 1:(xci.nMarkers*3)  %loop through all data columns

                [r,lags] = xcorr(xcm.data(k:wl+k,j),xci.data(k:wl+k,j),maxlag,'coeff');
                allrs(j,:) = r;

            end
          
          
          [c, p] = max(nanmean(allrs)); %TODO: try also to choose only the best half of the allrs curves (the ones with highest peak values)
          
          lagConfidences(end+1) = c;
                
          syncCandidates(end+1) = lags(p);
          
        end
          
        %This seems to be quite good...
        data(i).syncpoint = median(syncCandidates(lagConfidences > 0.1));

    end




end