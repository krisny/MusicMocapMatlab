function minOut = mcarrayMin(d,mintype)
% 
% Mean across all mocap structs in an array.
% If the number of frames in the structs are different, the longer structs
% will be cut.
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%

if nargin < 2
    mintype = 'normal';
end
    
    dframes = min([d.nFrames]);
    dmarkers = max([d.nMarkers]);

    dl = length(d);

    q = nan(dframes,dmarkers*3,dl);

    for i = 1:dl

        q(1:dframes,1:1:d(i).nMarkers*3,i) = d(i).data(1:dframes,:);

    end

    
    minOut = d(1);
    minOut.nFrames = dframes;
    minOut.nMarkers = dmarkers;
    minOut.filename = 'arraymin';
    
    if strcmpi(mintype,'absmin')

        
        [~,y] = nanmin(abs(q),[],3);
    
        x= (1:dframes*dmarkers*3)' + (y(:)-1)*(dframes*dmarkers*3);
        
        minOut.data = reshape(q(x(:)),dframes,dmarkers*3);
        
    else
        
        minOut.data = nanmin(q,[],3);
        
    end


end

