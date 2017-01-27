function mctimeautocorr(d)

colormap('default')
bluefromgray = gray;
bluefromgray(:,1)=0;
bluefromgray(:,2)=0;
m=[flipud(bluefromgray);hot];

mytempo=120; %calculated using mirtoolbox

[per, ac, ~, lag] = mcwindow(@mcperiod,d,300/mytempo,0.05);

tider = [11,53,60,74,88,102,117];

clf

colormap(m)

for i=1:d.nMarkers
subplot(d.nMarkers/2,2,i)
aci = (squeeze(ac(:,:,i)));
scalefactor = max(max(max(abs(ac))));
image(60/mytempo*(1:size(per,1)),lag,128*(0.5+aci/scalefactor));
%image(512*(0.5+aci/scalefactor));
set(gca,'YTick',0:30/mytempo:1.9,'YTickLabel',{'0','1/8','1/4','3/8','1/2','5/8','3/4','7/8','1/1'})
ylim([lag(6),lag(length(lag))])
title(['Marker ' d.markerName{i}])
ylabel('Periode')
xlabel('Tid i sekunder')

for j=1:length(tider)
    line([tider(j) tider(j)],[0 2],'color','w','LineStyle',':')
end

end

end


%7 11 21 41 51 81 101