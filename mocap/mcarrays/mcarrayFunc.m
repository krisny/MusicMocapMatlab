function dout = mcarrayFunc(d,f,varargin)

% Wrapper for applying native Mocap Toolbox functions to mcarrays
% dout = mcarrayFunc(data,f,varargin)
% 
% f: string, name of the function, followed by all arguments to the mocap
% toolbox function. 
%
% Examples:
%
% d2 = mcarrayFunc(data,'mcnorm')
% d2 = mcarrayFunc(data,'mctrim',2,7)
% d2 = mcarrayFunc(data,'mctrim',25,700,'frame')
% d2 = mcarrayFunc(data,'mctimeder',2)
% d2 = mcarrayFunc(data,'mcmean')
%
% By Kristian Nymoen, RITMO/University of Oslo, 2019
%


fh = str2func(f);

x = varargin;


if nargout(fh)
    switch class(fh(d(1),x{:}))
        case 'struct'
            for i = 1:length(d)
                dout(i) = fh(d(i),x{:});
            end
        otherwise
            for i = 1:length(d)
                dout{i} = fh(d(i),x{:});
            end
    end
else
    for i = 1:length(d)
        fh(d(i),x{:});
    end
end     
    
end