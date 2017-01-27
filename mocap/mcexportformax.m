function mcexportformax(d)



if isfield(d,'type') && max([strcmp(d.type, 'MoCap data') strcmp(d.type, 'norm data')])
    
    
    datafile1 = fopen(['./' d.filename '.txt'], 'w');
    
    %fprintf(datafile1, '%i ', [tmp(n,:)]); 
    fprintf(datafile1, 'Freq %i nMarkers %i dim %i %i timederOrder %i', [d.freq d.nMarkers size(d.data) d.timederOrder]); 
    fprintf(datafile1, '\n'); 

    d.dim = size(d.data,2)/d.nMarkers;
    
    for n = 1:size(d.data,1)
        
        for m = 1:d.nMarkers
            fprintf(datafile1, '/toMax/%s ', [d.markerName{m}]); 
            fprintf(datafile1, '%f ', [d.data(n,((m-1)*d.dim+1):((m-1)*d.dim)+d.dim)]); 
            fprintf(datafile1, '\b, '); 
        end
        fprintf(datafile1, '\b\b\b\n'); 
    end

    fclose(datafile1);

    clear n datafile1 tmp m 
    
    
    
 %below: things that work but not store as osc   
 if 0 
    datafile1 = fopen(['./' d.filename '.txt'], 'w');
    
    %fprintf(datafile1, '%i ', [tmp(n,:)]); 
    fprintf(datafile1, 'Freq %i nMarkers %i dim %i %i timederOrder %i', [d.freq d.nMarkers size(d.data) d.timederOrder]); 
    fprintf(datafile1, '\n'); 

    for n = 1:size(d.data,1)
        fprintf(datafile1, '%f ', [d.data(n,:)]); 
        fprintf(datafile1, '\n'); 
    end

    fclose(datafile1);

    clear n datafile1 tmp    
  end
else
    disp('The first input argument should be a variable with MoCap data structure.');
end


