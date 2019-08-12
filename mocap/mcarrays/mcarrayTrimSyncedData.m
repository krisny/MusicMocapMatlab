function data = mcarrayTrimSyncedData(data)

% data = mcarrayTrimSyncedData(data)
%
% input:  struct array of mocap structs, with a .syncpoint field. 
%         The syncpoint is given in frames, and may be obtained from the
%         mcarraySyncWindowed function
%
% output: struct array of mocap structs. Aligned according to the
%         syncpoint field, and cut by the latest start and earliest 
%         stop time across all recordings.
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%

    
    
    dnfold = [data.nFrames];
    padding = 500;


    for i = 1:length(data)

        offset(i) = round(data(i).syncpoint)+padding;
    
        data(i).data = [nan(offset(i)  , data(i).nMarkers*3); data(i).data; nan(offset(i)  , data(i).nMarkers*3);];
        data(i).nFrames =  data(i).nFrames + offset(i)*2; 
    
    end

    
    startTrim = max([offset])+1;
        
    endTrim = min([data.nFrames])-min(offset);


    for i = 1:length(data)

        data(i) = mctrim(data(i),startTrim,endTrim,'frame');
        data(i).syncpoint = 0;
    
    end



end