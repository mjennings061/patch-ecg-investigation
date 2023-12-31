function [ rWaveLocations ] = qrsDetector(waveformData, sampleFrequency)
% QRS Detection
% Based on Hamilton and Tompkins Method
%returns rWaveLocations which is the sample number of each r-wave detected.
%waveformData should be an vector of data (at least 10 dseconds in
%duration)
%sampleFruency is in Hz


Fs = sampleFrequency;
[stage4Out, stage1Out] = ECGPreprocessor(waveformData, sampleFrequency);
%rWaveLocations = 0;

%%%%%%%%%%%%Decission Logic Stage*************************
 TH = 0.189;                                                        %set coeffieicnt used in calculation of threshold.  Based on what is published

%set up buffers and variables to hold values during calculations of
%thresholds etc

QRS_Peak_Buffer = [];                                               %create a QRS peak buffer
for counter = 1:Fs:8*Fs                                             %fill buffer with which has max from 1 second blocks of the first 8 seconds of signal
    QRS_Peak_Buffer = [QRS_Peak_Buffer max(stage4Out(counter:counter+Fs-1))];
end
Noise_Peak_Buffer = zeros(1,8);                                      %set first 8 samples in noise buffer to 0
RR_interval_Buffer = ones(1,8);                                      %set first 8 samples in RR interval buffer to 1


%first build peak bufer using turning point methods in function
%"peakDetector"
bigPeakStore  = peakDetector( stage4Out );


%%%  now that all peaks are detected go through the peak array and classify
%%%  them as QRS peaks or noise
% mean(QRS_Peak_Buffer);
qrsThreshold = median(Noise_Peak_Buffer) + TH * (median(QRS_Peak_Buffer)- median(Noise_Peak_Buffer));
peakCat = 0;
qrsThresholdArray = [];
allFidVals = [];
oldDerivitave = 0;
oldFidLoc = 0;
actualFidLoc = 0;
z=1;                                                                    %so start at the start point in our data
previousPeakLoc = 0;    
i = 0;
j = 0;
k = 0;
searchBackFlag = 0;

while z <= length(bigPeakStore(:,1))                                    %go through our big peak store
    if (bigPeakStore(z,2) > qrsThreshold) && (((bigPeakStore(z,4) - previousPeakLoc)/Fs)>0.2) && (bigPeakStore(z,4) > round(Fs*.225))               %so if the peak vlau of the current peak is > qrs threshold and it has been 200ms since the last qrs also make sure there is no peak too close to the start of the record as this will cause a failure     

        %%we have established that this is actually a qrs complex and now
        %%we go off and now find the r-wave
        % we know it is a QRS as the recorded peak value of the turning
        % point is greater than the current QRS threshold
        % and the end of this peak and the previous peak is more than 0.2
        % seconds since the end of th elast peak and that we have seen at
        % least enough samples that will allow us to count back at least
        % .225 seconds later on
        %%%%%% 
        
        
        %In the below if statment we set up parameters thst relate to how
        %quickly the end of this peak appears.  If the peak occours long
        %before the end of the windowed signal then this may be a qrs and t
        %wave combined.  In that case we adjust the window that we search
        %to detect the r wave in the actual signal. In this case we also
        %set the end of the peak to be a new value (i.e. 170 ms after the
        %peak) and we ignore the end that we detected in the peak detect
        %section above.
        if (bigPeakStore(z,4) - bigPeakStore(z,3)) > (0.17*Fs)
            windowStart = 0.250;
            windowEnd = 0.150;
            averagedPeakEnd = bigPeakStore(z,3) + round(0.17*Fs);  
        else
            windowStart = 0.225;
            windowEnd = 0.125;
            averagedPeakEnd = bigPeakStore(z,4);
        end
        
stage1Out;%= waveformData;
        [~, fidLoc] = max(abs(stage1Out(averagedPeakEnd-round(Fs*windowStart):averagedPeakEnd-round(Fs*windowEnd))));       %now search the window from .225 to .150 seconds before the 50% of peakQRS max point 
        actualFidLoc = averagedPeakEnd-round(Fs*windowStart)+fidLoc-1;                                                      %now find the actual fid location by putting on the correct number of samples pror
        newDerivitave = max(diff(stage1Out(actualFidLoc-round(Fs*0.02):actualFidLoc+round(Fs*0.02))));                      %find the max of the derivative of the filtered signal for the 20mSec before the peak was detected        

        %the following section only proceeds to log the r wave detection if
        %there has been 360ms since the last detection or, if the new
        %deriviative value is more than half the value of the previous
        %derivative
        
        if  (((bigPeakStore(z,4) - previousPeakLoc)/Fs)>0.36) || (newDerivitave > 0.5*oldDerivitave)        %check to make sure the new derivative is greater than 0.5 of the previous. if it is not this is not a qrs (might be a t)
            allFidVals = [allFidVals; actualFidLoc];                                                        %store the fid location 
            peakCat = [peakCat 1];                                                                          %mark a flag to categorise this as a QRS peak
            QRS_Peak_Buffer = [bigPeakStore(z,2) QRS_Peak_Buffer];                                          %put the value in the start of our QRS bufer
            previousPeakLoc = bigPeakStore(z,4);                                                            %flag where our previous peak was located
            oldDerivitave = newDerivitave;                                                                  %set our derivitae to be the old derivitave
            currentRR = (actualFidLoc - oldFidLoc)/Fs;                                                      %now that we have a peak work out the rr interval
            RR_interval_Buffer = [currentRR RR_interval_Buffer];                                            %buffer the rr interval
            oldFidLoc =  actualFidLoc;                                                                      %keep the old fid location for calulating the rr later            
            searchBackFlag = 0;
        else
            peakCat = [peakCat 0];                                                                          %if the ocnditions have not been met record this as a noise peak and not a qrs
            Noise_Peak_Buffer = [bigPeakStore(z,2) Noise_Peak_Buffer];  
        end
        %before moving on we want to update the threshold and move on to
        %the next row in our big peak buffer
        
        if searchBackFlag == 0
            qrsThreshold = median(Noise_Peak_Buffer(1:8)) + (TH * (median(QRS_Peak_Buffer(1:8))- median(Noise_Peak_Buffer(1:8))));     
        end
        
        z = z+1;
        i = i+1;
        
    elseif (((bigPeakStore(z,4) - previousPeakLoc)/Fs) > RR_interval_Buffer(1)*1.5) && (previousPeakLoc + round(Fs/(1/0.36)) < bigPeakStore(z,4))

        %here we check and see if we have reached a point where we have not
        %had a QRS peak end for more than 150% of the rr interval
        

        %if so go back to the end of the last qrs peak and find where that is in
        %our big peak store. and move forward until we find a peak 360mS after this
        tempIndex = find(bigPeakStore(1:z,4) >= (previousPeakLoc + round(Fs/(1/0.36))));       
        %find the max peak that exists between that peak and now
        ((bigPeakStore(z,4) - previousPeakLoc)/Fs);
        (previousPeakLoc + round(Fs/(1/0.36)));
        [r,e] = max(bigPeakStore(tempIndex,2));
        tempPeakValueLoc = tempIndex(1)+e-1;
        if bigPeakStore(tempPeakValueLoc,2) > (0.5 * qrsThreshold)
        %if that peak exceeds 50% of the threshold then move the current
        %big peak pointer to that point and half the threshold, the peak
        %will then be picked up on the next cycle of the above code
        %we dont change the qrs peak bffer here
            qrsThreshold = 0.5 *qrsThreshold;
            z = tempPeakValueLoc;
            searchBackFlag = 1;
        else
            %if there is no peak excieeding 0.5 of the threshold then we do
            %nothing other than assume this is not a QRS peak and we
            %categorise the noise peak, get a new threshold on move the z
            %pointer on one.
            peakCat = [peakCat 0];
            Noise_Peak_Buffer = [bigPeakStore(z,2) Noise_Peak_Buffer]; 
            qrsThreshold = median(Noise_Peak_Buffer(1:8)) + (TH * (median(QRS_Peak_Buffer(1:8))- median(Noise_Peak_Buffer(1:8))));     
            z=z+1;
            j = j+1;
        end
        
    else
        %if we are in here this is not a QRS peak and we can assume that it
        %has not been more than 150% of rr since last peak.  we just move
        %on and calculate a new threshold and register this as a noise peak
        peakCat = [peakCat 0];
        Noise_Peak_Buffer = [bigPeakStore(z,2) Noise_Peak_Buffer]; 
        qrsThreshold = median(Noise_Peak_Buffer(1:8)) + (TH * (median(QRS_Peak_Buffer(1:8))- median(Noise_Peak_Buffer(1:8))));     
        z=z+1;
        k = k+1;
    end
    

end

rWaveLocations =  allFidVals;
bigPeaks = bigPeakStore;%[qrsThresholdArray];%peakCat';%[bigPeakStore];


end

function [ processedSignal, filteredECG ] = ECGPreprocessor(ECGdata, Fs)

lowerCutFreq = 3;
upperCutFreq = 18;
bandPassOrder = 2;
windowDuration = 0.160;

%%%%%%%%%Preprocessing Stage%%%%%%%%%%%%%%
%Stage 1
%Bandpass filter  
[b,a] = butter(bandPassOrder, [lowerCutFreq/(Fs/2) upperCutFreq/(Fs/2)], 'bandpass');
stage1Out = filtfilt(b,a,ECGdata);
filteredECG = stage1Out;

%Stage 2
%take first derivative using matlab function for differentiation 
stage2Out = diff(stage1Out);

%Stage 3 
%square signal
%check this as Hamilton may have proposed not to square here which is as
%opposed to what Pan and T did.
stage3Out = stage2Out.^2;

%Stage 4
%take a moving average of the signal.  Currently set to use a 160ms window
%please note moving average is calcualted looking back to previous samples
%as per original implmentation 

windowSampleSize = round(Fs*windowDuration);
stage4Out = zeros(length(stage3Out),1);
for z = windowSampleSize+1:length(stage3Out)
    stage4Out(z) = (1/windowSampleSize)*sum(stage3Out(z-windowSampleSize:z));
end

processedSignal =stage4Out;

end



function [ bigPeakArray ] = peakDetector( windowedSignal )


bigPeakStore = [];
%build peak array using turning point methods
z = 1;

while z < length(windowedSignal)
    if windowedSignal(z) > 0                                        %if we find any value greater than zero then we assume it is a peak
        zz = z;                                                     %if this is a peak we store the curent value of our index
        peakOnset = zz;                                             %also set a location for peak onset flag
        peakFlag = 1;                                               %now flag that we are in a peak region
        tempMaxPeakValue = windowedSignal(zz);                      %also set the first value of the peak to be the max value
        peakLocation = zz;                                          %and store the peak location
        while peakFlag == 1 && zz < length(windowedSignal)               %now move forward from this point (assuming we are in a peak) and that we haven't reached the end of the signal
            if windowedSignal(zz+1) >= tempMaxPeakValue                  %so check if the next value is greater than the current peak
                tempMaxPeakValue = windowedSignal(zz+1);                 %if it is set the temp peak vaue to the next value
                peakLocation = zz+1;                                %also mark the location of the new peak
                zz=zz+1;                                            %move our pointer on by one
            elseif windowedSignal(zz+1) <= (0.5*tempMaxPeakValue)        %now alternatively, if we have reached a value that is o.5 of the recorded peak
                peakOffset = zz+1;                                  %peak offset is set to th enext value (the o.5 peak max location)
                peakFlag = 0;                                       %indicate that we are out of a peak region
                zz=zz+1;                                            %we move our pointer on by one and soter the peak info (next line
                bigPeakStore = [bigPeakStore; [peakOnset tempMaxPeakValue peakLocation peakOffset]];
            else
                zz=zz+1;                                            %if we encouter somethng in the peak that isnt greater than the max or less than  half                                         
            end      
        end
        z = zz+1;                                                   %now set z to be the end of the last peak and then on the next line sotre allt he peak paramters 
    else
        z = z+1;
    end
end

bigPeakArray = bigPeakStore;

end