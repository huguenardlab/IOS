% code for capturing movies of IOS
if ~exist('tlCamera','var')
    [tlCameraSDK,tlCamera]=Prepare_camera;
end
% when done, be sure to close camera with this line:
% close_camera(tlCamera,tlCameraSDK)
clear pixels pixelsa;
% CAMERA PARAMETERS:
exposureTime=30; % exposure time for each frame in ms
frameRate=30; % in Hz
% Set exposure time and gain of the camera.
tlCamera.ExposureTime_us = exposureTime*1000; %convert from ms to us
tlCamera.FrameRateControlValue_fps=frameRate;
tlCamera.Disarm;
tlCamera.MaximumNumberOfFramesToQueue = 1; % just get the most recent frame and throw the rest away %numberOfFramesToAcquire;
tlCamera.Gain = 0; % GV modif 02/06
tlCamera.ROIAndBin.BinX=2; % JH modify 02/06/24 to speed imaging
tlCamera.ROIAndBin.BinY=2;
imageHeight = tlCamera.ImageHeight_pixels;
imageWidth = tlCamera.ImageWidth_pixels;
imagePixels=imageHeight*imageWidth;

figure(1);
axis image;
% Start software triggered image acquisition
disp('Starting software triggered image acquisition.');
% Set the number of frames per software trigger and start trigger
% acquisition
tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
% now go through and make the movies
tlCamera.Arm;
iloop=1;
tlCamera.IssueSoftwareTrigger;
% lineLength variable records number of characters in status message
lineLength=fprintf('Image frame number:%5d; %2d images queued',0,0);
notDone=1;
while notDone
    % Wait for image buffer to be filled to prevent sending too many
    % software triggers.
    while (tlCamera.NumberOfQueuedFrames == 0)
        pause(0.01);
    end
    % Get the pending image frame.
    imageFrame = tlCamera.GetPendingFrameOrNull;
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
        pixels(iloop)=thisFrame(320,256); % record one sample pixel intensity for this frame
        pixelsa(iloop)=mean(mean(thisFrame)); % record mean intensity for this frame
        iloop=iloop+1;
    end
    fprintf(repmat('\b',1,lineLength));
    fprintf('Image frame number:%5d; %2d images queued',imageFrame.FrameNumber,tlCamera.NumberOfQueuedFrames);
end

tlCamera.Disarm;
figure;plot(pixelsa);hold on; plot(pixels);title(sprintf("SDPixel= %.1f,SDFrame = %.1f, gain=%d",std(single(pixels)),std(single(pixelsa)),tlCamera.Gain));
ylim([500 700]);
