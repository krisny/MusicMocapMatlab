function meanOut = mcarrayMean(d)
% 
% Mean across all mocap structs in an array.
% If the number of frames in the structs are different, the longer structs
% will be cut.
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%

    
    dframes = min([d.nFrames]);
    dmarkers = max([d.nMarkers]);

    dl = length(d);

    q = nan(dframes,dmarkers*3,dl);

    for i = 1:dl

        q(1:dframes,1:1:d(i).nMarkers*3,i) = d(i).data(1:dframes,:);

    end

    meanOut = d(1);
    meanOut.nFrames = dframes;
    meanOut.nMarkers = dmarkers;
    meanOut.filename = 'arraymean';
    meanOut.data = nanmean(q,3);


end

