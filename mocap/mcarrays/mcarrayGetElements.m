function entries = mcarrayGetEntry(d,queryType,query,varargin)
% 
% Selects those entries of an array that satisfy a condition
%
%
% entries = mcarrayGetEntry(data,queryType,query,'name','value')
%
%
% Examples:
%  mcarrayGetEntry(data,'namecontains','Subject1')
%  mcarrayGetEntry(data,'namecontains','Subject1','outputtype','entries')
%  mcarrayGetEntry(data,'namecontains','subJECT1','ignoreCase',0)
%  mcarrayGetEntry(data,'nFrames',1000)
%  mcarrayGetEntry(data,'type','norm data')
%  mcarrayGetEntry(data,'type','segm data','outputtype','entries')
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%
%


outputType = 'data';
ignoreCase = 1;


for k=1:2:length(varargin)
    if strcmpi(varargin{k}, 'outputType')
        outputType=varargin{k+1};
    elseif strcmpi(varargin{k}, 'ignoreCase')
        ignoreCase=varargin{k+1};
    else        
        str=sprintf('Input argument %s unknown.', varargin{k});
        disp([10, str, 10])
    end
end




switch lower(queryType)
    case 'namecontains'
        [~,files,~]=cellfun(@fileparts,{d.filename},'Uni',0);
       
        if strcmpi(outputType,'indecies')
            entries = find(contains(files,query,'IgnoreCase',ignoreCase));
        else
            entries = d(contains(files,query,'IgnoreCase',ignoreCase));
        end
       

    otherwise
        if isfield(d,queryType)
            if isnumeric(query)
                if strcmpi(outputType,'indecies')
                    entries = find([d.(queryType)] == query);
                else
                    entries = d([d.(queryType)] == query);
                end
            else
                if strcmpi(outputType,'indecies')
                    entries = find(strcmpi({d.(queryType)}, query));
                else
                    entries = d(strcmpi({d.(queryType)}, query));
                end
            end
        else
            error('The second argument must be either ''namecontains'' or a field in the input struct')
        end
end

if isempty(entries)
    entries = [];
end




end

