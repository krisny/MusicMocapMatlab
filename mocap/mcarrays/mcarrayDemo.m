

% mcarray is a set of extension functions to the Mocap Toolbox for Matlab
% The toolbox is developed at the University of Jyväskylä, and can be found at:
% https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox
%
% The mcarray extensions allow batch processing of mocap recordings and
% statistical calculations across recordings.
%
% Developed by Kristian Nymoen, RITMO, University of Oslo, 2019
%
%
% This file takes you through the basic mcarray functions:
% 
% - mcarrayRead
% - mcarrayFunc
% - mcarrayGetEntry
% - mcarrayMax
% - mcarrayMean
% - mcarrayMedian
% - mcarrayMin
% - mcarrayPlotTimeseries
% - mcarrayStd
% - mcarrayMerge
% - mcarrayAnimate
% - mcarraySyncWindowed
% - mcarrayTrimSyncedData
%


%% Reading array data

% Grab your data with the mcarrayRead function, giving it
% the relative or absolute path to the mocap files as input. E.g.:
% filepath = '../../where/myData/is/'

data = mcarrayRead('./demoData');


%% 


%% Plotting
%  Like mcplottimeseries, the mcarrayPlotTimeseries requires that markernumber is specified. 

mcarrayPlotTimeseries(data,1) % marker #1 is the root marker, if 'dim' attribute is not set, then dimension 1 is plotted


%%
%using the same arguments as in mcplottimeseries, we may call for plotting
%markes 1 and 12, dimensions 1,2 and 3 like this:
mcarrayPlotTimeseries(data,[1 12],'dim',1:3) 

%%
% We'll return to plotting later, but let's close the plots before continuing
close all

%% Array functions
% You may run most of the native Mocap Toolbox functions on the data array
% using the mcarrayFunc function. Input the data array and the name of the 
% mocap toolbox function as arguments to the mcarrayFunc function.

%%
% First example, let's fill the gaps in all our recordings:

data = mcarrayFunc(data,'mcfillgaps');

%%
% some smoothing may be in order:

data = mcarrayFunc(data,'mcsmoothen');

%%
% let' also calculate the velocity and acceleration, and the norms of these

dataVel = mcarrayFunc(data,'mctimeder');
dataAcc = mcarrayFunc(data,'mctimeder',2); %note that arguments to the mctimeder are simply passed as additional arguments to mcarrayFunc
dataVelNorm = mcarrayFunc(dataVel,'mcnorm');
dataAccNorm = mcarrayFunc(dataAcc,'mcnorm');

%%
% and plot the velocity norm for root and feet markers in the same plot
mcarrayPlotTimeseries(dataVelNorm,[1 5,9],'plotopt','comb') 
%set xlim to only a short section of the recording:
xlim([20 22])



%% Syncing 
% The data we have been looking at so far are four recordings of jumping
% movement to a steady metronome. They are not synchronised to each other. 
% When the recordings contain the same type of movement to a steady
% metronome, we may attempt to synchronise them using the function
% mcarraySyncWindowed. The function uses cross correlation over a sliding
% window. 
%
% The default settings: data = mcarraySyncWindowed(data);
% would synchronise all other recordings to the first data entry in the 
% array, through a cross correlation across all dimensions of all markers. 
% The median best correlation is used as a syncpoint.
%

data = mcarraySyncWindowed(data);

% Note that the mcarraySyncWindowed function does not alter the data at
% all, but adds a .syncpoint field to the mocap data struct.

% Other settngs allows synchronising to a different data entry, to the mean
% in the array, changing window size, hop, maxlag, syncing by derivative,
% and selecting only a subset of markers. 
% See help mcarraySyncWindowed for more info

%%
% In order to align our recordings, we need to align trim them according to
% the extracted syncpoints. For this we use the mcarrayTrimSyncedData
% function.

data = mcarrayTrimSyncedData(datas);

% Let's also plot the third dimension of a selection of markers for our synced data:
mcarrayPlotTimeseries([datast],[1 5 9 12 16 20],'dim',3,'names',1);


%%
% The plots show us that the section between 13 and 28 seconds contain
% jumping motion for all participants, so let's trim our dataset further:

data = mcarrayFunc(data,'mctrim',13,28);

% We should also calculate new velocity and acceleration data:

dataVel = mcarrayFunc(data,'mctimeder');
dataAcc = mcarrayFunc(data,'mctimeder',2); %note that arguments to the mctimeder are simply passed as additional arguments to mcarrayFunc
dataVelNorm = mcarrayFunc(dataVel,'mcnorm');
dataAccNorm = mcarrayFunc(dataAcc,'mcnorm');


%% Returning to more plotting options
% plotting the new velocity norm data for root, head and left foot
mcarrayPlotTimeseries(dataVelNorm,[1 12 5],'plotopt','comb','names',1);

%%
% zooming in on a single jumping period
xlim([6 7])


%%
% As another plotting option, we can show the mean line for each marker
% across all recordings
close all
mcarrayPlotTimeseries(dataVelNorm,[1 12 5],'plotopt','comb','names',1,'plotMean',1);
xlim([6 7])

%%
% Or mean plus/minus standard deviation:
close all
mcarrayPlotTimeseries(dataVelNorm,[1 12 5],'plotopt','comb','names',1,'plotStd',1);
xlim([6 7])

%%
% The individual data lines may mess up the graph, these may be removed
% with the 'showlines' attribute:
close all
mcarrayPlotTimeseries(dataVelNorm,[1 12 5],'plotopt','comb','names',1,'plotStd',1,'showLines',0);
xlim([6 7])


%% Animating
% The mcarrayAnimate function displaces each of the entries in the data
% array and creates a video animation using mcanimate:

% First, let's make a shorter version of our data (just to save time)
shortdata = mcarrayFunc(data,'mctrim',7,10);

% Then, create an animation parameter structure, as usual for the mocap toolbox:
p = mcinitanimpar;
p.fps = 10;
p.conn = [1 2;2 3;3 4;4 5;1 6;6 7;7 8;8 9;1 10;10 11;11 12;11 13;13 14;14 15;15 16;11 17;17 18;18 19; 19 20];
p.msize= 3;
p.colors= 'kbwww';
p.output = 'test';
p.videoformat = 'mp4';

% And this is where the magic happens...
mcarrayAnimate(shortdata,p,0.6)


%%
% 3D animation may be generated by using the mc-variation functions at
% https://github.uio.no/krisny/mct-extras
%

p = mcinitanimpar('3D');
p.fps = 10;
p.conn = [1 2;2 3;3 4;4 5;1 6;6 7;7 8;8 9;1 10;10 11;11 12;11 13;13 14;14 15;15 16;11 17;17 18;18 19; 19 20];
p.msize= 3;
p.colors= 'wbwww';
p.output = 'test3D';
p.videoformat = 'mp4';
p.par3D.lightposition = [5000 20000 30000];
p.par3D.cameraposition = [15000 20000 4000];

mcarrayAnimate(shortdata,p,0.6)


%% Statistics
% As you have already seen in the plotting options, the mcarray tools can
% calculate certain simple statistics across the data entries

dataMean = mcarrayMean(shortdata);
dataMedian = mcarrayMedian(shortdata);
dataMin = mcarrayMin(shortdata);
dataMax = mcarrayMax(shortdata);

%% 
% Each of the calculated variables above is a single mc struct. Why don't
% we put them together and mcarrayAnimate them?

p.output = 'testStatistics';
mcarrayAnimate([dataMean dataMedian dataMin dataMax],p,0.6)


%% GetElements
% In cases where we need to look at only a subset our dataset, the 
% mcarrayGetElements may be helpful. 
%
% In our dataset, the file name contains the subject number, and we can
% extract a subset that matches a string in the file name like this:

dataS1 = mcarrayGetElements(data,'namecontains','01');

%% More advanced use of the mcarrayFunc
%
%
% The mcarrayFunc may also be used withd custom functions. To learn
% more about these, look at the functon help file: 
% >> help func
% Examples:

maxAccNormValue = mcarrayFunc(dataAccNorm,'@(d) max(max(d.data))');
sizeOfDataArray = mcarrayFunc(data,'@(d) size(d.data)');


