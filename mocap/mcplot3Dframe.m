function mcplot3Dframe(d, n, p, proj)
% Plots frames of motion capture data.  
% COMPATIBILITY NOTES (v. 1.5): Please use the function without the projection input argument, 
% but specify it in the animation structure instead.
%
% syntax
% par = mcplotframe(d, n);
% par = mcplotframe(d, n, p);
%
% input parameters
% d: MoCap data structure
% n: vector containing the numbers of the frames to be plotted
% p: animpar structure (optional)
% [depricated: proj: projection used: 0 = orthographic (default), 1 = perspective
%                    this flag is supposed to be set in the animation parameter stucture]
%
% output
% par: animpar structure used for plotting the frames (if color strings were used, they will converted to RGB triplets)
%
% examples
% par = mcplotframe(d, 1);
% mcplotframe(d, 500:10:600, par);
%
% comments
% If the animpar structure is not given, the function calls
% mcinitanimpar and sets the .limits field of the animpar structure
% automatically so that all the markers fit into all frames.
%
% see also
% mcanimate, mcinitanimpar
%
% 
% Based on mcplotframe in the Motion Capture Toolbox, 3D version developed
% by Kristian Nymoen, University of Oslo
% 
% ? Part of the Motion Capture Toolbox, Copyright ?2008,
% University of Jyvaskyla, Finland



par=[];

if isfield(d,'type') && strcmp(d.type, 'MoCap data') || isfield(d,'type') && strcmp(d.type, 'norm data') || isfield(d,'type') && strcmp(d.type, 'segm data')
else disp([10, 'The first input argument has to be a variable with valid mocap toolbox data structure.', 10]);
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    return;
end

if nargin>3 %BBADd0150303
    disp([10, 'Please note that, from MCT version 1.5, the perspective projection flag is set in the field "perspective" in the animation'])
    disp(['parameters. Please adapt your code accordingly.', 10])
    p.perspective=proj;
end


if nargin<3
    p = mcinitanimpar;
    if nargin==1 || ~isnumeric(n) %output fix [BB20111031]
        disp([10, 'Please set frame number(s) you want to plot.', 10]);
        [y,fs] = audioread('mcsound.wav');
        sound(y,fs);
        return
    end
end

for k=1:length(n)
    n1=n(k);
    if n1>d.nFrames
        w1=sprintf('Indicated frame(s) (%d) exceeds number of frames in data (%d).', max(n), d.nFrames);
        disp([10, w1, 10]);
        [y,fs] = audioread('mcsound.wav');
        sound(y,fs);
        return
    end
end

%for compatibility of parameter structure
if isfield(p, 'folder') %#BB20150303/20150717
    p.output=p.folder;
    p = rmfield(p,'folder');
    disp([10, 'Please note that, from MCT version 1.5, the parameter "output" is used instead of "folder" to specify the file/folder name.'])
    disp(['Please adapt your code accordingly.', 10])
end

if isfield(p,'type') && strcmp(p.type, 'animpar')
else disp([10, 'The third input argument has to be an animation parameter structure.', 10]);
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    return;
end 

if p.animate && p.getparams==0 %BBADd0150303
    currdir = cd; % store current directory
    if p.createframes==1
        mkdir(p.output) %BB20150507 in case frames (series of png files) are needed, and not video file
        cd(p.output) 
    else %VideoWriter is used with video file output %BB20150507
        fn=p.output; %BB_NEW_20140212 for VideoWriter
        if strcmp(p.videoformat,'avi')
            movObj = VideoWriter(fn); %set file name
        else
            movObj = VideoWriter(fn,'MPEG-4'); %BBADd0150717
        end
        movObj.FrameRate = p.fps; %set frame rate before opening the video object
        open(movObj); %open the object
    end
end




par=p;

% color management (compatibility): convert old string color definition into num array [BB20111031]
if ischar(p.colors)
    colors=NaN(5,3);
    for k=1:5
        colors(k,:)=lookup_l(p.colors(k));
    end
else colors=p.colors; %Fix BB20141202 ? if colors were already in num array
end

bgcol=colors(1,:); % set background color

if isfield(p,'markercolors') && ~isempty(p.markercolors) %field and colors are set
    if ischar(p.markercolors) % but in string format
        mcol=NaN(d.nMarkers,3);
        for k=1:size(p.markercolors,2)
            mcol(k,:)=lookup_l(p.markercolors(k));%convert to num array
        end
    else 
        mcol=p.markercolors; %field and colors are set in num format already
    end
    if d.nMarkers > length(p.markercolors)
        for k=length(p.markercolors)+1:d.nMarkers
            mcol(k,:)=colors(2,:);
        end
    end
else %no field or/and empty
    mcol=repmat(colors(2,:),d.nMarkers,1);
end

if isfield(p,'conncolors') && ~isempty(p.conncolors) %field and colors are set
    if ischar(p.conncolors) % but in string format
        ccol=NaN(size(p.conn,1),3);
        for k=1:size(p.conncolors,2)
            ccol(k,:)=lookup_l(p.conncolors(k));%convert to num array
        end
    else 
        ccol=p.conncolors; %field and colors are set in num format already
    end
    if size(p.conn,1) > length(p.conncolors)
        for k=length(p.conncolors)+1:size(p.conn,1)
            ccol(k,:)=colors(3,:);
        end
    end
else %no field or/and empty
    ccol=repmat(colors(3,:),size(p.conn,1),1);
end

%BBFIX 20120404: mcmerge problems
if p.trl~=0
    tcol=repmat(colors(4,:),d.nMarkers,1);
    if isfield(p,'tracecolors') && ~isempty(p.tracecolors) %(field and) tracecolors are set
        if ischar(p.tracecolors) % but in string format
            for k=1:size(p.tracecolors,2)
                tcol(k,:)=lookup_l(p.tracecolors(k));%convert to num array
            end
        else
            tcol(1:size(p.tracecolors,1),:)=p.tracecolors; %tracecolors in num format already
        end
    end
    if isempty(p.trm)
        p.trm=1:d.nMarkers; %plot all traces if trm is empty
    else
        if length(p.trm)<d.nMarkers
            %sort p.tracecolors/tcol1 according to p.trm
            %increase p.trm to have same length as d.nMarkers and fill it up with NaNs
            tmp=repmat(colors(4,:),d.nMarkers,1);
            tmp1=nan(d.nMarkers,1)';
            for k=1:length(p.trm)
                tmp(p.trm(k),:)=tcol(k,:);
                tmp1(p.trm(k))=p.trm(k);
            end
            tcol=tmp;
            p.trm=tmp1;
        end
    end
else
    tcol=[];
    p.trm=[];
end
p.tracecolors=tcol;

if ~isempty(p.trm)%created problem when merging and translating...
    p.trm=p.trm(:);
end

%warning if trace length is empty, but marker trace vector is set 
if (p.trl==0 && ~isempty(p.trm)) || (p.trl==0 && ~isempty(p.tracecolors))
    disp([10, 'Warning: Please set trace length (trl) in your animation parameters in order to plot traces.', 10])
end

%if trace length is set, but vector indicating the markers to be traced is empty, all markers will be trace
if p.trl~=0 && isempty(p.trm)
    disp([10, 'Note: All markers traced.', 10])
    p.trm=1:d.nMarkers;
end


%BBFIX20120404: mcmerge problems
if p.showmnum==1
    ncol=repmat(colors(5,:),d.nMarkers,1);
    if isfield(p,'numbercolors') && ~isempty(p.numbercolors) %(field and) numbercolors are set
        if ischar(p.numbercolors) % but in string format
            for k=1:size(p.numbercolors,2)
                ncol(k,:)=lookup_l(p.numbercolors(k));%convert to num array
            end
        else
            ncol(1:size(p.numbercolors,1),:)=p.numbercolors; %numbercolors in num format already
        end
    end
    if isempty(p.numbers)
        p.numbers=1:d.nMarkers; %plot all markers if numbers is empty
    else
        if length(p.numbers)<d.nMarkers
            %sort p1.numbercolors/ncol1 according to p1.numbers
            %increase p1.numbers to have same length as d1.nMarkers and fill it up with NaNs
            tmp=repmat(colors(5,:),d.nMarkers,1);
            tmp1=nan(d.nMarkers,1)';
            for k=1:length(p.numbers)
                tmp(p.numbers(k),:)=ncol(k,:);
                tmp1(p.numbers(k))=p.numbers(k);
            end
            ncol=tmp;
            p.numbers=tmp1;
        end
    end
else
    ncol=[];
    p.numbers=[]; %fill up p1.numbers with nan
end


par.colors=colors;
par.markercolors=mcol;
par.conncolors=ccol;
par.tracecolors=tcol;
par.numbercolors=ncol;

par.numbers=p.numbers;


%fill up trace widths (twidth) to have same length as traced markers (trm) / nMarkers
if p.trl~=0
    if length(p.twidth)<d.nMarkers
        twidth=nan(d.nMarkers,1)';
        i=1;
        for k=1:length(p.trm)
            if isnan(p.trm(k))
            else
                twidth(k)=p.twidth(i);
                if i<length(p.twidth)
                    i=i+1;
                end
            end
        end
        p.twidth=twidth;
    end
end
    
%individual widths for connectors and traces
%fill up cwidth
if length(p.cwidth)<size(p.conn,1)
    cl=length(p.cwidth);
    cwidth=nan(size(p.conn,1),1);
    cwidth(1:length(p.cwidth))=p.cwidth;
    for k=cl+1:length(cwidth)
        if isnan(cwidth(k))
            cwidth(k)=cwidth(cl);%fill up with last given value
        end
    end
elseif length(p.cwidth)>size(p.conn,1)%if cwidth is longer than conn
    cwidth=p.cwidth(1:size(p.conn,1));
else
    cwidth=p.cwidth;
end
p.cwidth=cwidth;

% if isfield(p,'cwidth') && length(p.cwidth)>1
%     p.cwidth=p.cwidth;
% else
%     p.cwidth=repmat(p.cwidth(1),1,length(p.conn));
% end
if isfield(p,'twidth') && length(p.twidth)>1
    p.twidth=p.twidth;
else
    if ~isempty(p.trm)
        p.twidth=repmat(p.twidth(1),1,length(p.trm));
    end
end
 

% az and el doesn't work that well in the 3d animation. using
% 'CameraPosition' attribute of axes object instead
%az=p.az;
%el=p.el;


%d1 = mcrotate(d, az, [0 0 1]);
%d = mcrotate(d1, el, [1 0 0]); %%%%%%%
    
if p.perspective==0 % orthographic projection    
    x=d.data(n,1:3:end);
    y=d.data(n,2:3:end);
    z=d.data(n,3:3:end);
else % perspective projection
    disp('perspective option disabled in 3D version, using orthographic')
    x=d.data(n,1:3:end);
    y=d.data(n,2:3:end);
    z=d.data(n,3:3:end);
%     if ~isfield(p,'pers') % for backward compatibility, use default values
%         p.pers.c=[0 -4000 0];
%         p.pers.th=[0 0 0];
%         p.pers.e=[0 -2000 0];
%     end
% 
%     th=180*p.pers.th/pi;
%     rot1=[1 0 0; 0 cos(th(1)) -sin(th(1)); 0 sin(th(1)) cos(th(1))];
%     rot2=[cos(th(2)) 0 sin(th(2)); 0 1 0; -sin(th(2)) 0 cos(th(2))];
%     rot3=[cos(th(3)) -sin(th(3)) 0; sin(th(3)) cos(th(3)) 0; 0 0 1];
%     dd=zeros(size(d.data(n,:)));
%     %p.pers.e(2)=min(min(d.data(n,2:3:end)));
%     for k=1:d.nMarkers
%         dd(:,3*k+(-2:0))=(rot1*rot2*rot3*(d.data(n,3*k+(-2:0))'-repmat(p.pers.c',1,length(n))))';
%     end
%     % make closest marker to be on the projection plan
%     dd(:,2:3:end)=dd(:,2:3:end)-min(min(dd(:,2:3:end)))+p.pers.e(2)-p.pers.c(2); 
%     x=-(dd(:,1:3:end)-repmat(p.pers.e(1),length(n),d.nMarkers)).*(p.pers.e(2)./dd(:,2:3:end));
%     y=(p.pers.e(2)-p.pers.c(2))./dd(:,2:3:end); % used for marker size scaling
%     z=-(dd(:,3:3:end)-repmat(p.pers.e(3),length(n),d.nMarkers)).*(p.pers.e(2)./dd(:,2:3:end));    
end


    tmp=d.data(:,1:3:end);tmp=tmp(:); maxx=nanmax(tmp); minx=nanmin(tmp);
    tmp=d.data(:,2:3:end);tmp=tmp(:); maxy=nanmax(tmp); miny=nanmin(tmp);%miny=miny-abs(miny*1.2);
    tmp=d.data(:,3:3:end);tmp=tmp(:); maxz=nanmax(tmp); minz=nanmin(tmp);
    midx = (maxx+minx)/2;
    midy = (maxy+miny)/2;
    midz = (maxz+minz)/2;
    
    maxx = maxx+abs(maxx*0.05);
    maxy = maxy+abs(maxy*0.05);
    maxz = maxz+abs(maxz*0.05);
    minx = minx-abs(minx*0.05);
    miny = miny-abs(miny*0.05);
    minz = minz-abs(minz*0.05);
    
    maxxyz = max([maxx,maxy,maxz]);


%campos = ones(1,3).*[maxx,maxy,maxz]*1.5
campos = [maxx,maxy,maxz].*[14 20 4]; %camera position
lightPos = [maxx,maxy,maxz].*[1.5 20 4];

% %BBADd0150303: exit function here without doing the animation or plotting, 
% but setting the parameters, esp. the limits, to make videos with a reduced 
% set of markers that look exactly like the videos with all markers
if p.getparams==1
    par=p;
    par = orderfields(par, {'type','scrsize','limits','az','el','msize','colors','markercolors',...
    'conncolors','tracecolors','numbercolors','cwidth','twidth','conn','conn2','trm','trl',...
    'showmnum','numbers','showfnum','animate','fps','output','videoformat','createframes','getparams','perspective','pers'});
    return
end

fignr=1;

if p.animate %20150720 / HJ: in animate case, set figure and axes outside main loop
    %figure(fignr); 
    clf;
    set(gcf, 'WindowStyle','normal');
    %set(gcf,'Position',[50 50 p.scrsize(1) p.scrsize(2)]) ; % DVD: w=720 h=420
    %set(gcf, 'color', bgcol);
    %view(0,90);
    colormap([ones(64,1) zeros(64,1) zeros(64,1)]);
end

if isfield(p,'floorimage')
    if ~isempty(p.floorimage)
        floorimg = imread(p.floorimage);
        floorlevel = minz;
        floorimgsize = size(floorimg); floorimgsize(end)=[];
        floorscale = max([maxx-minx,maxy-miny,maxz-minz])/min(floorimgsize);
    end
end
if isfield(p,'wallimage')
   if ~isempty(p.wallimage)
        wallimg = imread(p.wallimage);
        wallimgsize = size(wallimg); wallimgsize(end)=[];
        wallscale = max([maxx-minx,maxy-miny,maxz-minz])/min(wallimgsize);
        wallimg = flip(wallimg ,1);
   end
end

for k=1:size(x,1) % main loop
    if  p.animate
        clf; 
        axes('position', [0 0 1 1], 'XLim', [minx maxxyz], 'YLim', [miny maxxyz],'ZLim', [minz maxxyz], 'CameraPosition',campos,'CameraTarget',[midx midy midz],'Projection','perspective','CameraUpVector',[0 0 1],'Color',p.colors(1));
        %view([0 90])
        daspect([1 1 1])
        hold on;
    else
        figure(fignr); 
        clf;
        set(gcf, 'WindowStyle','normal');
        %set(gcf,'Position',[50 50 p.scrsize(1) p.scrsize(2)]) ; % DVD: w=720 h=420        
        axes('position', [0 0 1 1], 'XLim', [minx maxxyz], 'YLim', [miny maxxyz],'ZLim', [minz maxxyz], 'CameraPosition',campos,'CameraTarget',[midx midy midz],'Projection','perspective','CameraUpVector',[0 0 1],'Color',p.colors(1));
        daspect([1 1 1])
        hold on
        %set(gcf, 'color', bgcol);
        view(10,150);
        %colormap([ones(64,1) zeros(64,1) zeros(64,1)]);
        fignr=fignr+1;
    end
    
    %plot some text to appear in background
%     text(maxxx-650, minzz+350, {YOUR TEXT'}, 'FontSize', 24, 'color', [.15 .15 .15]);


%plot floor

 if isfield(p,'drawfloor')
     if p.drawfloor == 1
        floortransform = hgtransform('Matrix',makehgtform('xrotate',0,'scale',floorscale*1,'translate',[minx miny floorlevel]/floorscale*1));
        image(floortransform,floorimg)
     end
 end
    %plot walls
 if isfield(p,'drawwallx')
     if p.drawwallx == 1
        xbackwalltransform = hgtransform('Matrix',makehgtform('scale',wallscale,'translate',[minx miny minz]/wallscale,'xrotate',pi/2));
        image(xbackwalltransform,wallimg)
     end
 end
  if isfield(p,'drawwally')
     if p.drawwally == 1
        ybackwalltransform = hgtransform('Matrix',makehgtform('scale',wallscale,'translate',[minx miny minz]/wallscale,'xrotate',pi/2,'yrotate',pi/2));
        image(ybackwalltransform,flip(wallimg ,2))
     end
  end
  
    

        
  %image([0,0,0],[1000 10 1000],img)
%end
%view(3)
    
         
    % plot marker-to-marker connections
    if ~isempty(p.conn)
        if 0
        for m=1:size(p.conn,1)
            %if x(k,p.conn(m,1))*x(k,p.conn(m,2))~=0
            if isfinite(x(k,p.conn(m,1))*x(k,p.conn(m,2)))
                plot3([x(k,p.conn(m,1)) x(k,p.conn(m,2))], [y(k,p.conn(m,1)) y(k,p.conn(m,2))], [z(k,p.conn(m,1)) z(k,p.conn(m,2))], '-','Color',ccol(m,:),'LineWidth', p.cwidth(m));
                %plot3([x(k,p.conn(m,1)) x(k,p.conn(m,2))], [z(k,p.conn(m,1)) z(k,p.conn(m,2))],[y(k,p.conn(m,1)) y(k,p.conn(m,2))], '-','LineWidth', p.cwidth(m));
            end
        end
        else
            for m=1:size(p.conn,1)
                
%                c1 = mapar.conn(i,1);
%                c2 = mapar.conn(i,2);

                r1 = [x(k,p.conn(m,1)) y(k,p.conn(m,1)) z(k,p.conn(m,1))];
                r2 = [x(k,p.conn(m,2)) y(k,p.conn(m,2)) z(k,p.conn(m,2))];
                [pcx,pcy,pcz] = cylinder2P(p.cwidth*0.002*maxxyz,20,r1,r2);
                tmpbone = surf(pcx,pcy,pcz);
                tmpbone.EdgeColor = 'none';
                tmpbone.FaceColor = p.colors(3);
            end
            light('Position',lightPos)
        end
        
    end
    grid on
    % plot midpoint-to-midpoint connections
    if ~isempty(p.conn2)
        disp('conn2 not implemented yet')
       
%         for m=1:size(p.conn2,1)
%             %if x(k,p.conn2(m,1))*x(k,p.conn2(m,2))*x(k,p.conn2(m,3))*x(k,p.conn2(m,4))~=0
%             if isfinite(x(k,p.conn2(m,1))*x(k,p.conn2(m,2))*x(k,p.conn2(m,3))*x(k,p.conn2(m,4)))
%                 tmpx1 = (x(k,p.conn2(m,1))+x(k,p.conn2(m,2)))/2;
%                 tmpx2 = (x(k,p.conn2(m,3))+x(k,p.conn2(m,4)))/2;
%                 tmpy1 = (z(k,p.conn2(m,1))+z(k,p.conn2(m,2)))/2;
%                 tmpy2 = (z(k,p.conn2(m,3))+z(k,p.conn2(m,4)))/2;
%                 plot([tmpx1 tmpx2], [tmpy1 tmpy2], '-','Color',ccol(m,:),'LineWidth', p.cwidth(m));
%             end
%         end

    end
    
    % plot traces if animation
    if p.animate && p.trl~=0
%         trml=sort(p.trm);%BB: sorting marker traces 
        trlen = round(p.fps * p.trl);
        start=max(1,k-trlen);
        ind = start-1+find(~isnan(x(start:k)));
        for m=1:length(p.trm)
            if isnan(p.trm(m)) %BBFIX 20120404 mcmerge adaption - NaN traces not plotted 
            else
                plot3(x(ind,m),y(ind,m),z(ind,m),'-','Color',tcol(m,:),'Linewidth',p.twidth(m));
            end
        end
    end
    
    
    % plot markers
    if 0
    for m=1:size(x,2)
        %if x(k,m)~=0 & ~isnan(x(k,m)) % if marker visible
        if isfinite(x(k,m)) % if marker visible
                plot3(x(k,m),y(k,m),z(k,m),'o','MarkerSize',p.msize(min(m,length(p.msize))),'MarkerEdgeColor',mcol(m,:),'MarkerFaceColor','k')
            if p.showmnum
                if isempty(p.numbers)
                    h=text(x(k,m),y(k,m),z(k,m)+maxxyz/80,num2str(m));
                    set(h,'FontSize',16);
                    set(h,'Color',ncol(m,:))
                else
                    
                     %if ismember(m, p.numbers) %FIX BB20120326: redefinining numbers plotting (mcmerge) 
                     if m<=length(p.numbers)
                         if isnan(p.numbers(m)) %NaN numbers not plotted
                         else
                             h=text(x(k,m),y(k,m),z(k,m)+maxxyz/50,num2str(p.numbers(m)));
                             set(h,'FontSize',16);
                             set(h,'Color',ncol(m,:))
                         end
                    end
                end
            end
%            axis off
        end


    end
    else
       
        [px,py,pz] = sphere(50);                % generate coordinates for a 50 x 50 sphere

        px=px*p.msize*0.002*maxxyz;
        py=py*p.msize*0.002*maxxyz;
        pz=pz*p.msize*0.002*maxxyz;
        
        for m=1:size(x,2)

            sEarth(m) = surface(px+x(k,m), py+y(k,m),flip(pz)+z(k,m));   
            sEarth(m).FaceColor = p.colors(2); 
            sEarth(m).EdgeColor = 'none';              % remove surface edge color
            %sEarth(i).CData = floorimg;                   % set color data 
            
            SP = shadowPoint([0 1 0],[0 miny 0],lightPos,[x(k,m) y(k,m) z(k,m)]);
            line(SP(1),SP(2),SP(3),'Marker','o','MarkerFaceColor','k','MarkerSize',4,'color','k')
            SP = shadowPoint([1 0 0],[minx 0 0],lightPos,[x(k,m) y(k,m) z(k,m)]);
            line(SP(1),SP(2),SP(3),'Marker','o','MarkerFaceColor','k','MarkerSize',4,'color','k')
            SP = shadowPoint([0 0 1],[0 0 minz],lightPos,[x(k,m) y(k,m) z(k,m)]);
            line(SP(1),SP(2),SP(3),'Marker','o','MarkerFaceColor','k','MarkerSize',4,'color','k')

            
        end

        
        
    end
    
    
    
    
    %axis off
    
    if p.showfnum
        h=text(minx+0.95*(maxx-minx), miny+0.05*(maxy-miny), minz+0.05*(maxz-minz), num2str(k),...
            'HorizontalAlignment','Right','FontSize',12,'FontWeight','bold');
        set(h,'Color',colors(5,:))
    end
    
%     %plot some copywrite text or anything else - BB20121102
%      text(maxxx-650, minzz+350, {'Jyv?skyl? Music & Motion Capture'}, 'FontSize', 24, 'color', [.9 .9 .9]);
%     text(maxxx-300, 0, {'Birgitta Burger', 'Jyv?skyl? Univ.', 'Finland'});
%     text(minxx+40, minzz+0.97*(maxzz-minzz), 'High Sub-Band 2 Flux', 'FontSize', 16, 'FontWeight', 'bold');
%     text(minxx+70, minzz+0.97*(maxzz-minzz)-75, {'high speed of head'}, 'FontSize', 12, 'FontWeight', 'bold');
    
    
    drawnow
    hold off
    if p.animate
        if p.createframes==1 
            fn=['frame', sprintf('%0.4d',k),'.png']; %old version: create frames
            imwrite(frame2im(getframe),fn,'png');
%             fn=['frame', sprintf('%0.4d',k),'.eps'];
%             saveas(gcf, fn, 'eps');
        else
            writeVideo(movObj,getframe(gcf)); %BB_NEW_20140212 for VideoWriter
        end
    end
end

if p.animate
    if p.createframes==0
        %close
        close(movObj);
    else %close(gcf);
    end
    cd(currdir);
end





return;

end


function colorar=lookup_l(colorstr)
if strcmp(colorstr, 'k')
    colorar=[0 0 0];
elseif strcmp(colorstr, 'w')
    colorar=[1 1 1];
elseif strcmp(colorstr, 'r')
    colorar=[1 0 0];
elseif strcmp(colorstr, 'g')
    colorar=[0 1 0];
elseif strcmp(colorstr, 'b')
    colorar=[0 0 1];
elseif strcmp(colorstr, 'y')
    colorar=[1 1 0];
elseif strcmp(colorstr, 'm')
    colorar=[1 0 1];
elseif strcmp(colorstr, 'c')
    colorar=[0 1 1];
end

return;

end

function SP = shadowPoint(planeNormalVec,pointOnPlane,p1,p2)

    SP = p1+(-dot(planeNormalVec,p1 - pointOnPlane) / dot(planeNormalVec,p2-p1)).*(p2-p1);

end




