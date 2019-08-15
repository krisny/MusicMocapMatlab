function data = mcarrayRead(datafolder)

% Read all the mocap data in a folder into a struct array
% mcarrayRead(datafolder)
% 
% default datafolder = './data/'
%

    if nargin == 0
        datafolder = ['.' filesep 'data' filesep]; %look in subfolder called 'data'
    end

    if ~strcmp(datafolder(end),filesep)
        %require a slash at the end of the folder path
        datafolder = [datafolder filesep];
    end

    files = dir([datafolder '*.tsv']);
    if isempty(files)
        files = dir([datafolder '*.mat']);
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
        disp(['Warning: No mocap files found in ' datafolder])
    end

    
    for i = 1:length(files)
        
        data(i) = mcread([datafolder files(i).name]);
        
    end
    

end
