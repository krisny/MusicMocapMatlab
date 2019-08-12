function data = mcarrayRead(datafolder)

% Read all the mocap data in a folder into a struct array
% mcarrayRead(datafolder)
% 
% default datafolder = ./data/
%

    if nargin == 0
        datafolder = './data/';
    end


    files = dir([datafolder '*.mat']);
    if isempty(files)
        files = dir([datafolder '*.tsv']);
    end
    if isempty(files)
        files = dir([datafolder '*.c3d']);
    end
    if isempty(files)
        files = dir([datafolder '*.bvh']);
    end
    if isempty(files)
        files = dir([datafolder '*.wii']);
    end

    if isempty(files)
        disp('no mocap files found in datafolder...')
    end

    
    for i = 1:length(files)
        
        data(i) = mcread([datafolder files(i).name]);
        
    end
    

end