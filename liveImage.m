% code for capturing movies of IOS
if ~exist('tlCamera','var')
    [tlCameraSDK,tlCamera]=Prepare_camera;
end
% when done, be sure to close camera with this line:
% close_camera(tlCamera,tlCameraSDK)
addpath 'Q:\matlabscripts\WhiskerStim'


clear pixels pixelsa;
% CAMERA PARAMETERS:
exposureTime=30; % exposure time for each frame in ms
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


%Set the FIFO frame buffer size. Default size is 1.
tlCamera.Disarm;
movieDuration=3.5; % seconds
numberOfFramesToAcquire = movieDuration*frameRate;
tlCamera.MaximumNumberOfFramesToQueue = 1; % just get the most recent frame and throw the rest away %numberOfFramesToAcquire;

numberOfMovies=20;
timeBetweenMovies=20; %seconds


% Check if the camera supports setting "Gain"
gainRange = tlCamera.GainRange;
if (gainRange.Maximum > 0)
    tlCamera.Gain = 3;
end
tlCamera.Gain = 0; % GV modif 02/06
tlCamera.ROIAndBin.BinX=2; % JH modify 02/06/24 to speed imaging
tlCamera.ROIAndBin.BinY=2;
imageHeight = tlCamera.ImageHeight_pixels;
imageWidth = tlCamera.ImageWidth_pixels;
imagePixels=imageHeight*imageWidth;

% FILES SETUP
secondsPerDay=60*60*24;
filePrefix=datestr(now,'mm-dd-yy');
fileSerialNumber=1;
dataFilePath='Q:\\Camille\IOS\';
if ~isfolder([dataFilePath filePrefix])
    mkdir([dataFilePath filePrefix]);
end
dataFilePath=[dataFilePath filePrefix '\'];



% Set the FIFO frame buffer size. Default size is 1.
%tlCamera.MaximumNumberOfFramesToQueue =  50; % this seems to be enough not to drop frames, which are now kept in results.frametimes
figure(1);
axis image;
% Start software triggered image acquisition
disp('Starting software triggered image acquisition.');
% Set the number of frames per software trigger and start trigger
% acquisition
tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
%hFig=figure (1);
 %       set(hFig,'WindowKeyPressFcn',@keyPressCallback);

% now go through and make the movies
    tlCamera.Arm;
    iloop=1;
    tlCamera.IssueSoftwareTrigger;
      lineLength=fprintf('Image frame number:%5d; %2d images queued',0,0);
    notDone=1;
    while notDone
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
        maxInt=40;
        

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
            thisFrame=flipud(reshape(uint16(imageData), [imageWidth,imageHeight]));
            maxIntensity=max(max(thisFrame));
            minIntensity=min(min(thisFrame));
            meanIntensity=mean(mean(thisFrame));
            sdIntensity=std(std(single(thisFrame)));
            sd=0;
            if iloop>2
                sd=std(single(pixels));
            end
                figure(1),imagesc(thisFrame), colormap(gray), colorbar ; %,caxis([0 1023])
                title(sprintf('Frame %3d: Mean/Min/Max Int.= %.1f/%d/%d, SD1px=%.1f SDfrm=%.1f',iloop,meanIntensity,minIntensity, maxIntensity,sd,sdIntensity));
                pixels(iloop)=thisFrame(320,256);
                pixelsa(iloop)=mean(mean(thisFrame));
                
                %pixels(iloop)=mean(mean(thisFrame));
                
                iloop=iloop+1;
                   % set(gcf,'WindowKeyPressFcn',@keyPressCallback);
if iloop >199
   % notDone=0; % uncomment this to stop imaging after 200 frames if
   % desired
end
               % imagesc(reshape(uint16(imageData), [imageWidth,imageHeight])');
                %colorbar;
                %caxis([0 maxIntensity]);
            end
            fprintf(repmat('\b',1,lineLength));
            fprintf('Image frame number:%5d; %2d images queued',imageFrame.FrameNumber,tlCamera.NumberOfQueuedFrames);
          end
    tlCamera.Disarm;

figure;plot(pixelsa);hold on; plot(pixels);title(sprintf("SDPixel= %.1f,SDFrame = %.1f, gain=%d",std(single(pixels)),std(single(pixelsa)),tlCamera.Gain));
ylim([500 700]);
   