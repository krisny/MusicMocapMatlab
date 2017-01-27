function mcspectrogram(d)

colormap('default')
% bluefromgray = gray;
% bluefromgray(:,1)=0;
% bluefromgray(:,2)=0;
% m=[flipud(bluefromgray);hot];

%mytempo=120; %calculated using mirtoolbox

[f, sp] = mcwindow(@mcspectrum,d,1000,0.05);

%tider = [11,53,60,74,88,102,117];

clf

% colormap(m)

for i=1:d.nMarkers
    subplot(d.nMarkers,1,i)
    spi = (squeeze(sp(:,:,i)))';
    scalefactor = 1;%max(max(max(abs(sp))));
    image(64*(spi/scalefactor));
    %image(64*(spi/scalefactor));
    set(gca,'YTick',0:5:100,'YTickLabel',100:-5:0)
    ylim([85,100])
    title(['Marker ' d.markerName{i}])
    %ylabel('Periode')
    %xlabel('Tid i sekunder')
    
%for j=1:length(tider)
%    line([tider(j) tider(j)],[0 100],'color','w','LineStyle',':')
%end

end

end


%7 11 21 41 51 81 101