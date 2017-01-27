function corrs = knwindowedxcorr(a,b,duration,ws,hop)
% Creates a plot of the correlation coefficients of cross correlated, windowed streams of dat
%
% syntax
% knwindowedxcorr(a,b);				
% knwindowedxcorr(a,b,duration);                Specifies the duration of the two data streams (in seconds) - for axis labeling purpose only
% knwindowedxcorr(a,b,duration,windowsize);     Specifies window size (in samples). Default 1/100 of duration
% knwindowedxcorr(a,b,duration,windowsize,hop);	Specifies hop size (between 0 and 1). Default 0.1
%
% Developed by Kristian Nymoen, University of Oslo, 2014



if nargin < 3
    duration = 1;
end

if nargin < 5
    hop=0.1;
end


colormap('default')
bluefromgray = gray;
bluefromgray(:,1)=0;
bluefromgray(:,2)=0;
m=[flipud(bluefromgray);hot];


if length(a)>length(b)
    a=resample(a,length(b),length(a));
elseif length(b)>length(a)
    b=resample(b,length(a),length(b));
end

a=a-mean(a);
b=b-mean(b);



if nargin < 4
    ws = floor(length(a)/100); %if no window size given: set ws = 1/100 of length
end

corrs = [];

for i = 1:ws*hop:(length(a)-ws)


%    [c,lags] = xcorr(a(i:(i+ws-1)),b(i:(i+ws-1)),'unbiased');
%    [c,lags] = xcorr(a(i:(i+ws-1)),b(i:(i+ws-1)),'biased');
    [c,lags] = xcorr(a(i:(i+ws-1)),b(i:(i+ws-1)),'coeff');
%    [c,lags] = xcorr(a(i:(i+ws-1)),b(i:(i+ws-1)),'biased');
    corrs = [corrs,c];
   

end


clf

subplot(4,1,1:3)

colormap(m)

scalefactor = 1;%max(max(abs(corrs)));

image(0:duration/size(corrs,1):duration,duration/length(a)*lags,128*(0.5+corrs/scalefactor));

ylabel('Lag in seconds')

title(sprintf('Cross correlations for window size: %.2g seconds, Hop: %.2g', duration/length(a)*ws, hop))

subplot(4,1,4)



[maxs,ind]=max(mean(corrs,2));

%Plotte korrelasjon ved 0 lag
steps=0:((duration)/(size(corrs,2)-1)):duration;

%Plotting 
plot(steps,corrs(ind,:));
%title(['Highest average coorelation is at ' num2str(duration/length(a)*lags(ind)) ' seconds:'])

title(sprintf('Highest average correlation is %.2g at lag %.3f seconds:', max(maxs),duration/length(a)*lags(ind)))

%Plotting at 0 lag:
%plot(steps,corrs(size(corrs,2)/2,:));

xlim([0 duration])

if duration == 1
    xlabel('Time (normalised)')
else
    xlabel('Time (sec)')
end

ylabel('correlation')




end


