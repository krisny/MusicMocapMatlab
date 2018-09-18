function dist = mcmarkerdist2(d, m1, m2)
% Calculates the frame-by-frame distance of a marker pair.
%
% syntax
% dist = mcmarkerdist(d, m1, m2);
%
% input parameters
% d: MoCap data structure
% m1, m2: marker numbers
%
% output
% dist: column vector
%
% examples
% dist = mcmarkerdist(d, 1, 5);
%
% (Updates by Kristian Nymoen, University of Oslo. Work in progress:
% - m2 may be an array of two markers, giving the distance between marker 
%   m1 and a line spanned by two other markers
% - m2 may be an array of three markers, giving the distance between 
%   marker m1 and a plane spanned by the three other markers
%
% ? Part of the Motion Capture Toolbox, Copyright ?2008, 
% University of Jyvaskyla, Finland

if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    if length(m2) == 1
        %calc distance between two markers
        if m1>d.nMarkers || m2>d.nMarkers || m1<1 || m2<1 
            disp('Marker numbers out of range.'); 
            return; 
        end
        c1 = 3*m1+(-2:0); c2 = 3*m2+(-2:0); 
        dist = sqrt(sum((d.data(:,c1)-d.data(:,c2)).^2,2));
    elseif length(m2) == 2
        %calc distance between marker and line segment
        if m1>d.nMarkers || max(m2>d.nMarkers) || m1<1 || max(m2<1)
            disp('Marker numbers out of range.'); 
            return; 
        end
        c1 = 3*m1+(-2:0); c2 = 3*m2(1)+(-2:0); c3 = 3*m2(2)+(-2:0);
        %I'm sure there is a more elegant way of doing this, but at least it works:
        for i = 1:d.nFrames
            dist(i) = norm(cross([d.data(i,c2);d.data(i,c3)]',[d.data(i,c2);d.data(i,c1)]'))/norm([d.data(i,c2);d.data(i,c3)]');
        end
    elseif length(m2) == 3
        %calc distance between marker and plane
        if m1>d.nMarkers || max(m2>d.nMarkers) || m1<1 || max(m2<1)
            disp('Marker numbers out of range.'); 
            return;
        end
        c1 = 3*m1+(-2:0); c2 = 3*m2(1)+(-2:0); c3 = 3*m2(2)+(-2:0); c4 = 3*m2(3)+(-2:0);
        %another ugly for loop solution... 
        for i = 1:d.nFrames
            c = cross(d.data(i,c2)-d.data(i,c3),d.data(i,c2)-d.data(i,c4));
            n = c/sqrt(sum(c.^2));
            dist(i) = n*(d.data(i,c1)-d.data(i,c2))';
        end
    else
        disp('The third argument must be either one, two (line), or three (plane) markers');
    end
else
    disp('The first input argument should be a variable with MoCap data structure.');
end

