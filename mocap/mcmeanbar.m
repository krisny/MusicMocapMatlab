
function md = mcmeanbar(d,barlength,smoothFrontalMarkers)
% 
% Frame-by-frame average across successive fixed periods
%
% mcmeanbar(d,barlength)
% mcmeanbar(d,barlength,smoothFrontalMarkers)
%
% d: mocap data structure
% barlength: repetition length (in frames)
% smoothFrontalMarkers: optional
%       uses a smoothed version of mc2frontal to ensure the subject
%       maintaines the general same direction. 
%       A cell array of marker numbers denoting the forward facing
%       direction is required. e.g. {[2 15];[7 19]};
%
% meandata = mcmeanbar(data,400,{[2 15];[7 19]});
%


if nargin > 2

    tm2j = mcinitm2jpar;
    tm2j.nMarkers = 2;
    tm2j.markerNum = smoothFrontalMarkers;
    tm2j.markerName = {'tmpL';'tmpR'};
    tmp = mcsmoothen(mcm2j(d,tm2j),[2, 0.01]);
    tmp = mcconcatenate(d,1:22,tmp,1:2);
    d = mc2frontal(tmp,23,24,'frame');
    d = mcconcatenate(d,1:22);
end

clear tmp

md = d;
md.nFrames = barlength;
nbars = (d.nFrames-1)/barlength;

for j = 1:nbars
    
    tmp(:,:,j) = d.data(((j-1)*barlength+1):(j*barlength),:);
    
end

md.data = nanmean(tmp,3);


end
