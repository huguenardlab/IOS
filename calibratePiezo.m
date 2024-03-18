% code for capturing movies of IOS
if ~exist('tlCamera','var')
    [tlCameraSDK,tlCamera]=Prepare_camera;
end
% when done, be sure to close camera with this line:
% close_camera(tlCamera,tlCameraSDK)
addpath 'Q:\matlabscripts\WhiskerStim'

prompt = {'Enter Mouse ID:','Enter Whisker:','Enter Piezo #','Enter User'};
dlgtitle = 'Whisker Stim Run';
fieldsize = [1 45; 1 45;1 45; 1 45];
definput = {'0000','D1','1','Gaby'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);


% CAMERA PARAMETERS:
exposureTime=20; % exposure time for each frame in ms
frameRate=30; % in Hz
% Set exposure time and gain of the camera.
tlCamera.ExposureTime_us = exposureTime*1000; %convert from ms to us
tlCamera.FrameRateControlValue_fps=frameRate;
% Set the FIFO frame buffer size. Default size is 1.
% if this is large enough, then we should not have any dropped frames.
% we can confirm this later in terms of the frame serial number saved along
% with each frame.  If the frame numbers are the same as the serial numbers
% for that frame, then we did not drop any.  Similary, the frame times
% should be separated by the reciprocal of the frame rate.  If there are
% any dropped frames, then frame times that are different than the next
% expected time as calculated by frame rate, then that is also evidence of
% dropped frames.  In tests where tlCamera.MaximumNumberOfFramesToQueue is 
% small, then both measures, frametimes and framenumbers DO indeed show 
% this discrepancy= numberOfFramesToAcquire;

tlCamera.Disarm;
%Set the FIFO frame buffer size. Default size is 1.
movieDuration=2.5; % seconds
numberOfFramesToAcquire = movieDuration*frameRate;
tlCamera.MaximumNumberOfFramesToQueue = numberOfFramesToAcquire;
numberOfMovies=5;
timeBetweenMovies=10; %seconds
imageHeight = tlCamera.ImageHeight_pixels;
imageWidth = tlCamera.ImageWidth_pixels;
imagePixels=imageHeight*imageWidth;

movieBuffer=zeros(numberOfFramesToAcquire,imagePixels);

% Check if the camera supports setting "Gain"
gainRange = tlCamera.GainRange;
if (gainRange.Maximum > 0)
    tlCamera.Gain = 0;
end


% FILES SETUP
secondsPerDay=60*60*24;
filePrefix=datestr(now,'mm-dd-yy');
fileSerialNumber=1;
dataFilePath='q:\a_LabResources\WhiskerStim\PiezoCalib\';
baseFilePath=dataFilePath;
if ~isfolder([dataFilePath filePrefix])
    mkdir([dataFilePath filePrefix]);
end
dataFilePath=[dataFilePath filePrefix '\'];


%WHISKERSTIM PARAMETERS
stim_cal=80; %80 um/V  This is the output scaling to use for the Polytec Vibrometer.  Also set tracking to slow and periodically hit displacement reset

desiredmovment = 300; %this is um micromovment desire
probe1calibration=30; %72 um/volts
stimScaler=desiredmovment/probe1calibration;
stimScaler=.125;
stim_increment=1;
total_pulses = 22;              % Number of pulses
pulse_frequency=10;
pulse_period=1/pulse_frequency;
pulse_duration = 0.050;         % Pulse duration in seconds (58 ms)
rise_time = 0.004;              % Rise time in seconds (4 ms)
delay_duration = 0.2;           % Delay before pulses start in seconds
sampleFrequency = 10000;            % Sampling rate in Hz
sweepLengthSeconds=movieDuration;
%stimScaler=1; % this is the calibration value for whisker stimular device 
 % the value in volts to apply to the Piezo driver to get the correct
 % amplitude, which should be ~187 um 5 mm from skin, which equals 2
 % degrees total deflection. our devices have a max of 7 volts drive
 
stimBuffer=WhiskerTemplate(total_pulses, pulse_frequency, pulse_duration, delay_duration, rise_time,sampleFrequency,sweepLengthSeconds,stim_increment)*stimScaler;% with 5 seconds of duration per sweep, and sampling at 10khz
totalSamples=sweepLengthSeconds*sampleFrequency;
stimDelay=sweepLengthSeconds/10;
dq=daq("ni"); % this opens any NI devices
%CameraTriggerOut=addoutput(dq,'Dev1','port1/line0','digital')
%write(dq,[1]); %This turns the output on (value =1, which is 5 V).
fprintf("max stim V %.2f",max(stimBuffer));
NIDev='Dev1';
dq.Rate=sampleFrequency;
CameraFrames=addinput(dq,NIDev,'ai0','Voltage');
StimOut=addoutput(dq,'Dev1','ao0','Voltage');
figure(1);
subplot(2,3,1);
h=plot(0);
title('Live data');
ylabel('Position (\mum)');
xlabel('Time (s)');
yl=[-100;400];
ylim(yl);
subplot(2,3,4);
h1=plot(0);
ylabel('Whisker Stim');
xlabel('Time (s)');

yl=[-.1;5.1];
ylim(yl);

% Set the FIFO frame buffer size. Default size is 1.
%tlCamera.MaximumNumberOfFramesToQueue =  50; % this seems to be enough not to drop frames, which are now kept in results.frametimes
subplot(2,3,3);
% Start software triggered image acquisition
disp('Starting software triggered image acquisition.');
% Set the number of frames per software trigger and start trigger
% acquisition
tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
todaysRuns=dir([dataFilePath '*']);
currentRun=0;
for f=1:size(todaysRuns,1)
    if strfind(todaysRuns(f).name,'Run')
        temp=str2num(extractAfter(todaysRuns(f).name,'Run'));
        if temp>currentRun
            currentRun=temp;
        end
    end
end
currentRun=currentRun+1;
runName=['Run' int2str(currentRun)];
if ~isfolder([dataFilePath runName])
    mkdir([dataFilePath runName]);
end
dataFilePath=[dataFilePath  '\' runName '\'];
clear allData;

% now go through and make the movies
for m=1:numberOfMovies
    todaysfiles=dir([dataFilePath filePrefix '*struct.mat']);
    currentserial=0;
    for f=1:size(todaysfiles,1)
        temp=str2num(todaysfiles(f).name(10:13));
        if temp>currentserial
            currentserial=temp;
        end
    end
    currentserial=currentserial+1;
    fileSerialNumber=currentserial;
    thisFileName=sprintf('%s%s-%04d',dataFilePath,filePrefix,fileSerialNumber);

    starttime=now;
    nextTime=starttime;
    nextTime=nextTime+timeBetweenMovies/secondsPerDay;

    flush(dq);
    preload(dq,stimBuffer');
    tlCamera.Arm;
    fprintf('Acquiring Data for %.1f seconds\n',sweepLengthSeconds);
    start(dq); %,"Duration", seconds(SweepLengthSeconds));
    
    tlCamera.IssueSoftwareTrigger;
    data2=[];
      lineLength=fprintf('Image frame number:%5d; %2d images queued',0,0);
    
    for iloop = 1:numberOfFramesToAcquire
        
        if  dq.Running
            data1=read(dq,'all','OutputFormat', 'Matrix');
            prevcollected=size(data2,1);
            data2=[data2; data1];
            collected=size(data1,1);
            if collected>=100
                set(h,'XData',((1:collected)+prevcollected)/sampleFrequency,'YData',data1(:,1)*stim_cal);
                set(h1,'XData',((1:collected)+prevcollected)/sampleFrequency,'YData',stimBuffer((1:collected)+prevcollected));
            end
            pause(.001);
        else
            try
                data1=read(dq,'all','OutputFormat', 'Matrix');
                data2=[data2; data1];
            catch
            end
        end

        % Wait for image buffer to be filled to prevent sending too many
        % software triggers.
        while (tlCamera.NumberOfQueuedFrames == 0)
            pause(0.01);
        end
        % If data processing in Matlab falls behind camera image
        % acquisition, the FIFO image frame buffer could be filled up,
        % which would result in missed frames.
%        if (tlCamera.NumberOfQueuedFrames > 1)
%            ll2=fprintf('%2d frames queued.     \n',tlCamera.NumberOfQueuedFrames);
%            lineLength=fprintf('Image frame number:%5d',0);
%        end
        % Get the pending image frame.
        playingMovie=0;
        imageFrame = tlCamera.GetPendingFrameOrNull; % this gets the actual time in ns of the frames, so we can make sure they
        % are evenly spaced at the frame rate.  each should be separated by
        % 33 ms, and they seem to be.  it looks like we should set this
        % value very high 
        if ~isempty(imageFrame)
            imageTimeStamp=imageFrame.TimeStampRelative_ns_OrNull.Value;
            frameNumber(iloop)=imageFrame.FrameNumber;
            frameTime(iloop)=imageTimeStamp;
            % For color images, the image data is in BGR format.
            imageData = imageFrame.ImageData.ImageData_monoOrBGR;
            if playingMovie
                figure(1);
                hh=subplot(2,3,2);
                cla(hh);
                colormap(gray);
                imagesc(reshape(uint16(imageData), [imageWidth, imageHeight])');
            end
            movieBuffer(iloop,:)=uint16(imageData);
            fprintf(repmat('\b',1,lineLength));
            fprintf('Image frame number:%5d; %2d images queued',imageFrame.FrameNumber,tlCamera.NumberOfQueuedFrames);
          end
    end
       try
                data1=read(dq,'all','OutputFormat', 'Matrix');
                data2=[data2; data1];
            catch
       end
    data=data2';
    allData(:,m)=data;
    figure(2); plot(mean(allData,2)*5120);
    figure(1);
    tlCamera.Disarm;
    frameTime2=double(frameTime-frameTime(1))/1e9;
    movieBuffer = reshape(movieBuffer, [numberOfFramesToAcquire,imageWidth, imageHeight]);
    if ~playingMovie
    hh=subplot(2,3,2);
    cla(hh);
    colormap(gray);
    fprintf('\n');
    lineLength= fprintf('Playing movie, Frame %5d, time %.2f s',i,i/frameRate);
    for i=1:3:numberOfFramesToAcquire
        imagesc(squeeze(movieBuffer(i,:,:))');
        fprintf(repmat('\b',1,lineLength));
        fprintf('Playing movie, Frame %5d, time %.2f s',i,i/frameRate);
        pause(1/frameRate);
        drawnow;
        %colorbar;
        title(['Frame' int2str(i)]);
    end
    fprintf('\n');
    end
    
    baselineEndsSeconds=1;
    lastBaselineFrame=baselineEndsSeconds*frameRate;
    baselineFrameAve=squeeze(mean(movieBuffer(1:lastBaselineFrame,:,:)));
    stimFrameAve=squeeze(mean(movieBuffer(lastBaselineFrame+1:end,:,:)));
    activityAve=stimFrameAve./baselineFrameAve;
    subplot(2,3,5)
    activityAveStack(m,:,:)=activityAve;
    imagesc(activityAve');
    title('Delta from Background');
    colormap(gray);
    colorbar('Location','southoutside');
    if m>1
        subplot(2,3,6);
        imagesc(squeeze(mean(activityAveStack))');
        title(sprintf("Delta averaged across %d movies",m));
        colormap(gray);
        colorbar('Location','southoutside');
    end
    subplot(2,3,3);
    hold off;
    plot((1:size(data,2))/sampleFrequency,stimBuffer(1:1:size(data,2)));
    hold on;
    plot((1:size(data,2))/sampleFrequency,data);
    legend({'WhiskerStim';'CameraFrames'});
    title(sprintf('Frames/stims - frame rate %.2f, %d dropped',tlCamera.GetMeasuredFrameRate,max(frameNumber)-numberOfFramesToAcquire));
    xlabel('Time (s)');
    lineLength= fprintf('%6.1f seconds until next loop',(nextTime-now)*secondsPerDay);
    save([dataFilePath runName '-diffImage'],"activityAveStack");
    results.pulse_frequency=pulse_frequency;
    results.pulse_duration = pulse_duration;
    results.rise_time=rise_time;
    results.delay_duration = delay_duration;
    results.data=data;
    results.frameTimes=frameTime2; % each frame time should be separated from the next by the 1/frame rate.
    % If it is bigger, then it should be a multiple, which will indicate number of skipped frames
    results.frameNumbers=frameNumber;  %if there are no dropped frames, each frame recorded should be the same as the one collected
   
    results.sampleFrequency=sampleFrequency;
    results.timeBetweenMovies=timeBetweenMovies;
    results.numloops=numberOfMovies;
    results.mouseID=answer{1};
    results.Whisker=answer{2};
    results.sweepLengthSeconds=sweepLengthSeconds ;
    save([thisFileName '-struct'],'results');
    while (nextTime>now)
        fprintf(repmat('\b',1,lineLength));

        fprintf('%6.1f seconds until next loop',(nextTime-now)*secondsPerDay);
        pause(.1);
    end

end
figure(3);
meanData=mean(allData')*stim_cal;
baseline=mean(meanData(1:10));
meanData=meanData-baseline;
timePoints=(1:size(data,2))/sampleFrequency;
plot(timePoints,meanData);
stims=find(stimBuffer>0);  % any point in the stim buffer that has a positive value, i.e. during a pulse
stimEnds=find(diff(stims)>1); % find the ends of each pulse

%stimEnds=[stimEnds stims(end)];
endSamples=([stims(stimEnds) stims(end)])-rise_time*sampleFrequency; 
stimEndTimes=timePoints(endSamples);
for i=1:size(stimEndTimes,2)
    movement(i)=mean(meanData(endSamples(i)-10:endSamples(i)));
    stim(i)=stimScaler*i;
end
hold on;
plot(timePoints(endSamples),movement,'o');
ylabel('position (\mum)');
xlabel('Time (s)');
figure(4);

plot(stim,movement,'o');
ylabel('position (\mum)');
xlabel('Drive (V)');
slope=stim'\movement'   ;
piezoNumber=str2num(answer{3});
user=answer{4};

title(sprintf('Piezo # %d. Drive: Slope %.1f \\mum/V',piezoNumber,slope));
intercept=movement(1)-slope*stim(1);
hold on;
plot(stim,stim*slope);
piezoCalibFile=[baseFilePath '/PiezoCalib.csv'];
if ~isfile(piezoCalibFile)
    f=fopen(piezoCalibFile,'w');
    fprintf(f,'PiezoNumber,CalibSlope(um/V),user,date\n');
    fclose(f);
end
dateStr=datestr(now);
    f=fopen(piezoCalibFile,'a');
    temp=fprintf(f,'%d,%.1f,%s,%s\n',piezoNumber,slope,user,dateStr);
    fclose(f);
save([dataFilePath '\Probe-' num2str(piezoNumber) '-calibration.mat'],'stim','movement');
