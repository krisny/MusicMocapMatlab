function d3 = mcappend (d1, d2)
% appending one mocap struct after another in time

if min(strcmp(d1.markerName,d2.markerName))
    d3 = mcinitstruct;
    d3.data = [d1.data;d2.data];
    d3.nFrames = d1.nFrames+d2.nFrames;
    d3.markerName = d1.markerName;
    d3.nMarkers = d1.nMarkers;
    d3.freq = d1.freq;
    d3.residual = [];
else
    error('Marker names of input data structs do not match!')
end
    
    



end