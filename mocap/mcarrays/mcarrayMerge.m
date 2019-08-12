function [allData, allPar] = mcarrayMerge(d,p,f,type,n)
% 
% Merge and translate array of mocap data structs
% mcarrayMerge(d)
% mcarrayMerge(d,p)
% mcarrayMerge(d,p,f)
% mcarrayMerge(d,p,f,type)
% mcarrayMerge(d,p,f,type,n)
% 
% INPUT:
% d: array of mocap data structs
% p: animation param struct (same for all entries in d)
% f: spatial separation between merged data entries
% type: translation type. 'circle' (default) positions the data entries
%       evenly in a circle, 'grid' aligns the entries in a grid with n
%       columns
% n: the number of grid columns
% 
% OUTPUT:
% allData: single mocap data struct (merged array)
% allPar: single animation param struct with bone connections for all data entries
%
% by Kristian Nymoen, University of Oslo
%


l = min([d.nFrames]);
if length(unique([d.nFrames])) > 1
    for i = 1:length(d)
        d(i) = mctrim(d(i),1,l,'frame');
    end
end

dmax = max(max([d.data]));
dmin = min(min([d.data]));
drange = dmax-dmin;

dl = length(d);

if nargin < 2
    p = [];
end

if nargin < 3
    f = 1;
end

if nargin < 4
    type = 'circle';
end

if nargin < 5
    n = 6;
end


for i = 1:dl
    
    if strcmpi(type,'circle')
        a = mctranslate(d(i),[cos(i*2*pi/dl)*drange*f*dl/20 sin(i*2*pi/dl)*drange*f*dl/20 0]); 
    else
        
        a = mctranslate(d(i),[mod(i-1,n)*drange*f double(idivide(int16(i-1),n))*drange*f 0]); 
    end
    

    if i == 1
        allData = a;
        allPar = p;
    else
        if isempty(p)
            allData = mcmerge(allData,a); 
        else
            [allData,allPar] = mcmerge(allData,a,allPar,p);
        end
    end
    
end
    

end