function mcarrayAnimate(q,apar,spacing,gridsize)

%input a struct array of mocap recordings. optionally with animation
%parameter struct, and specifying the layout of all recordings in space.
%
% use:
% - mcarrayAnimate(q,apar,spacing,gridsize)
% - mcarrayAnimate(q,apar)
% - mcarrayAnimate(q,apar,0.6,4)
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%


if nargin < 2
    apar = mcinitanimpar;
end

if nargin < 3
    spacing = 1;
end

if nargin < 4
    gridsize = 10;
end
    

[a,apar2] = mcarrayMerge(q,apar,spacing,'grid',gridsize);


mcanimate(a,apar2)



end