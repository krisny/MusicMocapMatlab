function d3 = mcsubtract(d1,d2,fno)

%calculate the difference between two mocap structs 

    if nargin < 3
        fno = 0;
    end
    
    
    
    d3 = d1;
    
    if fno
        d3.data = d1.data - d2.data(fno,:);
    else
        d3.data = d1.data - d2.data;
    end


end