%Return a median ECG beat from a 15-beat window
%9-lead beats (V1-V6, I-III)

function beatOutput  = averageBeat(ecgsigAll, plotFlag)
    MIN_BEAT_LENGTH = 500; %minimum length of the beats in samples
    MIN_WINDOW_LENGTH = 9000;   
    MAX_WINDOW_LENGTH = 25000;
    MIN_GAP = 500; %minimum number of samples between beats
    beatOutput = [];
    
    %% TEST - SPACIAL VELOCITY AND VCG
    order = [7 8 1 2 3 4 5 6];
    for i = 1:8
        ecg8(:,i) = ecgsigAll(:,order(i));
    end
    
    %% Filtering
    %band pass filter design
    fs=1000; %sampling freq (consider removing to outside function)
    h  = fdesign.bandpass('N,F3dB1,F3dB2', 2, 0.05, 150, fs); %2nd order 0.05-150Hz filter
    bpFilt = design(h, 'butter');
    %zero-padding before filtering
    N = length(ecgsigAll(:,1));
    N_leads = length(ecgsigAll(1,:));
    pad = zeros(round(N/2),N_leads);
    ecgsigAll = [pad;ecgsigAll;pad];
    ecgsigAllFilt = zeros(length(ecgsigAll),length(ecgsigAll(1,:)));
    for j = 1:length(ecgsigAll(1,:))    %for each lead (9 in total)
        %band pass filter 0.05-150Hz each lead
        ecgsigAllFilt(:,j) = filtfilt(bpFilt.sosMatrix, bpFilt.ScaleValues, ...
                                    ecgsigAll(:,j)); %zero-phase filtered ECG
    %         ecgsigAllFiltNotch(:,j) = filtfilt(notchFilt.Numerator, notchFilt.Denominator, ...
    %                                     ecgsigAllFilt(:,j)); %zero-phase notch filtered ECG
    end
    %remove zero padding post-filtering
    ecgsigAllFilt = ecgsigAllFilt(length(pad)+1:length(ecgsigAllFilt(:,1))-length(pad),:);

    %% Find peaks/R-waves using RMS beat
    ecgRMS = zeros(1,length(ecgsigAllFilt(:,1)));   %preallocation
    for sample = 1:length(ecgsigAllFilt(:,1)) %calculate rms of each sample
        ecgRMS(sample) = sqrt(mean(ecgsigAllFilt(sample,:).^2)); %RMS beat
    end
    rWaves = qrsDetector(ecgRMS, fs);     % use the pan-thompkins algorithm to detect R-waves

    %% Extract 15 consecutive beats by wavelet decomp
    %Select 15 beat window (first beat-gap between previous beat/2 : last beat + gap between next beat/2)
    % if there is a noise spike between any beat of the window, disregard and move to the next
    % window.
    % This is currently using d3 only, but this can be expanded to the higher frequencies
    ecgV2 = ecgsigAllFilt(:,2); %lead V2 is used for the noise calculations
    n_adjust = 15; %Change this to increase nSamples after r-wave
    for j = 1:length(rWaves)  %for all detected beats
        noiseFlag = 0;  %if this is 1, the window is rejected due to noise
        if(j<20)                %start on the 20th beat
            continue;
        end
        w_startBeat = rWaves(j-14); %get the first R-wave
        w_endBeat = rWaves(j);      %get the last R-wave (+15 beats later)
        if(j==15)   %the first beat start-point will be the first sample
            w_nStart = 1;
        else w_nStart = w_startBeat - round((w_startBeat-rWaves(j-15))/2) - n_adjust; %halfway between current beat and previous beat
        end
        if(j==length(rWaves))        %the last beat
            w_nEnd = length(ecgV2);    % set to last sample
        else w_nEnd = w_endBeat + round((rWaves(j+1)-w_endBeat)/2) + n_adjust; %halfway between last beat and next beat
        end
        w_ecg = ecgV2(w_nStart:w_nEnd); %ECG signal window (15 beats)

        %Find local peaks in window
%         [w_pks, w_locs] = findpeaks(abs(w_ecg).^2,'MinPeakHeight',0.2,'MinPeakDistance',MIN_GAP);
        %Fetch pks and locs for the window
        w_locs = rWaves(j-14:j)-w_nStart;
        if(length(w_locs) ~= 15) %if there's not 15 beats, ignore it
            clear w_ecg
            continue;
        end
        if(length(w_ecg) < MIN_WINDOW_LENGTH || length(w_ecg) > MAX_WINDOW_LENGTH) %if the window is too small or large, slide along
            clear w_ecg
            continue;
        end
        %Wavelet decompose window
        [cWin,l] = wavedec(w_ecg,5,'db4');
%         approxWin = appcoef(cWin,l,'db4');
        d5Win = wrcoef('d',cWin,l,'db4',5);
        d4Win = wrcoef('d',cWin,l,'db4',4);
        d3Win = wrcoef('d',cWin,l,'db4',3);
        d2Win = wrcoef('d',cWin,l,'db4',2);
        d1Win = wrcoef('d',cWin,l,'db4',1);

        %For each beat, check if noise on d3-d5 is above threshold (+/- 20 samples)
        N_DELAY = 40;
        for beat = 2:length(w_locs)  %for each of the 15 beats in the window w_ecg
            if(w_locs(beat)+N_DELAY > length(w_ecg)) %if its the last beat and it overflows
                maxBeat1_d3 = max(abs(d3Win(w_locs(beat)-N_DELAY:end))); %the maximum beat amplitude (current beat)
            else
                maxBeat1_d3 = max(abs(d3Win(w_locs(beat)-N_DELAY:w_locs(beat)+N_DELAY))); %the maximum beat amplitude (current beat)
            end

            if(w_locs(beat-1)-N_DELAY <= 0) %if its the first beat and it will overflow
                maxBeat2_d3 = max(abs(d3Win(1:w_locs(beat-1)+N_DELAY))); %the maximum beat amplitude (beat before)
            else
                maxBeat2_d3 = max(abs(d3Win(w_locs(beat-1)-N_DELAY:w_locs(beat-1)+N_DELAY))); %the maximum beat amplitude (beat before)
            end
            maxBeat = max([maxBeat1_d3 maxBeat2_d3]); %check which one has the higher amplitude
            if(~all(d3Win(w_locs(beat-1)+N_DELAY:w_locs(beat)-N_DELAY) < maxBeat)) %if the gap between beats has noise above maxBeat:
                noiseFlag = beat;                                          %set the noise flag and break
                break;
            end
        end

        %Accept and store window if its not noisy, continue to next step (else, repeat for next beat)
        if(noiseFlag > 0)
            noiseFlag = 0; %reset the noise flag
            continue;
        else
            %% Extract individual beats from window
            w_ecg9 = ecgsigAllFilt(w_nStart:w_nEnd, :); %9-lead beats (V1-V6, I-III)
            % Expand to 12-lead ECG
            %w_ecg12 = [w_ecg9 zeros(length(w_ecg9(:,1)),3)]; %12-leads (V1-V6,I-III, aVL,aVR,aVF)
            N_DELAY = 15; %Number of samples before/after the halfway point to adjust for
            beat_indv = cell(1,length(w_locs)); % Individual beats prellocation
            for beat = 1:length(w_locs) %for every beat in w_locs
                if(beat == 1) %if its the first beat and it will overflow
                    beat_indv{beat} = w_ecg9(1:w_locs(beat) ...
                        + round((w_locs(beat+1)-w_locs(beat))/2) + N_DELAY,:);
                elseif(beat == length(w_locs)) %if its the last beat
                    beat_indv{beat} = w_ecg9(w_locs(beat) ... 
                        - round((w_locs(beat)-w_locs(beat-1))/2) - N_DELAY:end,:);
                else
                    beat_indv{beat} = w_ecg9(w_locs(beat) ... 
                        - round((w_locs(beat)-w_locs(beat-1))/2) - N_DELAY : ...
                        w_locs(beat) + round((w_locs(beat+1)-w_locs(beat))/2) + N_DELAY,:);
                end
            end

            %% Align the beats based on v2
            %for all beats
            beat_lengths = zeros(1,length(beat_indv)); %preallocation
            for beat = 1:length(beat_indv)
                %find minimum length of all beats
                beat_lengths(beat) = length(beat_indv{beat}(:,1)); %length of beat 
            end

            beat_min = min(beat_lengths); %minimum beat length
            %for all beats   
            for beat = 1:length(beat_indv)
                %trim off the difference in length from the start of each beat
                diff = length(beat_indv{beat}(:,1)) - beat_min; %difference in beat lengths
                if(diff > 0)
                    beat_indv{beat}(1:diff,:) = []; %trim off the difference
                end
            end

            %get the first r-peak location from V2
            r_n = zeros(1,length(beat_indv)); %preallocation
            delta_r = zeros(1,length(beat_indv)); %preallocation
            [~,r_n(1)] = max(abs(beat_indv{1}(:,2))); %crude as hell. Works for now
            %for each beat after first beat (2:15)
            for beat = 2:length(beat_indv)
                %get peak sample no
                [~,r_n(beat)] = max(abs(beat_indv{beat}(:,2)));
                %circularly shift to align R 
                delta_r(beat) = r_n(1) - r_n(beat);
                beat_indv{beat} = circshift(beat_indv{beat}, delta_r(beat));
            end

            %trim the start and end of each beat by the maximum deviation between beats
            trim = round(max(delta_r));
            if( ((length(beat_indv{1}(:,1))-trim) < 1) || ...
                    (length(beat_indv{1}(:,1)) <= trim*2) ) %check if the beat size is too small
                disp(['WARNING: ignored due to small beat size'])
                clear w_ecg beat_indv
                continue;
            end
            for beat = 1:length(beat_indv)
                beat_indv{beat}(1:trim,:) = []; %trim off the difference
                beat_indv{beat}(end-trim:end,:) = []; %trim off the difference
            end

            %% Combine beats to an average beat for each lead
            beat_indv_mat = zeros(length(beat_indv{1}(:,1)), length(beat_indv{1}(1,:)), length(beat_indv)); %preallocation
            for beat = 1:length(beat_indv)
                beat_indv_mat(:,:,beat) = beat_indv{beat}; %convert from cell to 3D matrix
            end
            beatMean = zeros(length(beat_indv{1}(:,1)), length(beat_indv{1}(1,:))); %preallocation
            for lead = 1:length(beat_indv{1}(1,:))          %for each lead (1:9)
                for sample = 1:length(beat_indv{1}(:,1))    % and each sample (1:x)
    %                     beatRMS(sample,lead) = sqrt(mean(beat_indv_mat(sample,lead,:).^2)); %RMS beat
                     beatMean(sample,lead) = median(beat_indv_mat(sample,lead,:)); %calculate mean beat for each lead
                end
            end

            %Sanity check the beat length
            if(length(beatMean(:,1)) < 500)
                beatMean = [];
                disp(['WARNING: average beat ignored due to small beat size'])
                clear w_ecg beat_indv beat_indv_mat beatMean
                continue;
            end

            %% Return the average beat
            beatOutput = beatMean;

            %% Plot data 
            if(plotFlag == 1)
                %plot recording with 4-peaks for V2
                figure('Position',[0 50 1500 300])
                lead = 2; %lead V2=2 (V1-V6; I-III)
                plot(ecgsigAllFilt(:,lead))
                hold on;
                stem(rWaves,zeros(1,length(rWaves))); %r-peak locations
                grid on;
                xlabel('Sample No.');
                ylabel('Amplitude (mV)');
                title(['Lead ', num2str(lead), ' with R-peak locations']);
                legend('Signal','R locations');
                hold off;

                %plot chosen window - wavelet decomposition
                figure('Position',[0 50 800 750])
                subplot(6,1,1)
                plot(w_ecg)
                hold on
                plot(w_locs, 0, 'Marker', 'v','Color', 'g');
                hold off;
                title('Filtered ECG Window')
                subplot(6,1,2)
                plot(d5Win)
                title('Level 5 Detail Coefficients')
                subplot(6,1,3)
                plot(d4Win)
                title('Level 4 Detail Coefficients')
                subplot(6,1,4)
                plot(d3Win)
                title('Level 3 Detail Coefficients')
                subplot(6,1,5)
                plot(d2Win)
                title('Level 2 Detail Coefficients')
                subplot(6,1,6)
                plot(d1Win)
                title('Level 1 Detail Coefficients')
                sgtitle('Wavelet Decomposition of ECG Window (db4)')
                hold off;

                %plot overlayed beat
                figure()
                lead = 2; %lead V2=2 (V1-V6; I-III)
                for beat = 1:length(beat_indv)
                    plot(beat_indv{beat}(:,lead));
                    hold on;
                end
                grid on;
                xlabel('Sample number');
                ylabel('Amplitude (mV)');
%                 title(['Aligned beats of patient ', num2str(i), ', lead ', num2str(lead)]);
                title(['Aligned beats of patient, lead ', num2str(lead)]);
                hold off;

                %plot average beats 
                figure('Position',[800 50 800 750])
                subplot(3,3,1);
                plot(beatMean(:,7));
                grid on;
                title('I');
                subplot(3,3,2);
                plot(beatMean(:,1));
                grid on;
                title('V1');
                subplot(3,3,3);
                plot(beatMean(:,4));
                grid on;
                title('V4');
                subplot(3,3,4);
                plot(beatMean(:,8));
                grid on;
                title('II');
                subplot(3,3,5);
                plot(beatMean(:,2));
                grid on;
                title('V2');
                subplot(3,3,6);
                plot(beatMean(:,5));
                grid on;
                title('V5');
                subplot(3,3,7);
                plot(beatMean(:,9));
                grid on;
                title('III');
                subplot(3,3,8);
                plot(beatMean(:,3));
                grid on;
                title('V3');
                subplot(3,3,9);
                plot(beatMean(:,6));
                grid on;
                title('V6');
                grid on;
%                 sgtitle(['Average Beats for Patient ', num2str(i)])
                sgtitle('Average Beats for Patient')
%                 disp(['Patient ', num2str(i), ' outputs graphed']);
                disp('Outputs graphed');
            end
            break; %break the pks for loop and move to next patient
        end
    end
%     clear N N_leads pad ecg ecgsigAll ecgsigAllFilt ecgV2 w_ecg w_locs ...
%         w_endBeat w_nEnd w_nStart beat maxBeat2_d3 maxBeat1_d3 pad ...
%         w_ecg9 w_ecg12 approxWin cWin d1Win d2Win d3Win d4Win d5Win ...
%         beat_indv beat_lengths beat_min beat noiseFlag w_startBeat ...
%         beat_indv_mat r_n delta_r trim beatMean sample lead locs diff ...
%         N_DELAY n_adjust maxBeat l j pks;
    if(isempty(beatOutput)) %check if there is a valid output
        disp('WARNING: no average beat detected - Signal may be too noisy or too short');
%         w_nStart = 0;
%         w_nEnd = 0;
    end
end

function vcgLeads = calculateVCG(twelveLead_8)
%derive VCG 12 lead ECG
%8 independant leads of 12 lead required (I II V1 V2 V3 V4 V5 V6)
%VCG returnded as [XYZ]
%in both cases the rows are samples columns are leads
%Currently uses Guldenring
vcgXform = [0.5169 -0.2406 -0.0715;
-0.0722 0.6344 -0.1962; 
-0.0753 0.1707 -0.4987; 
0.0162 -0.0833 -0.0319; 
0.0384 0.1182  -0.2362; 
0.0545 0.0237 -0.0507; 
0.1384 -0.1649 -0.2007; 
0.4606 0.2100 0.4122];
 
vcgLeads = twelveLead_8 * vcgXform;
 
end
 
function spatialV = spatialVelocity(threeLeads,sampleFreq)
%calcualte spatial velocity signal
%returnded calues is square root of sum of swuares of first derivative of
%three leads
 
spatialV = sampleFreq*sqrt((diff(threeLeads(:,1)).^2)+(diff(threeLeads(:,2)).^2)+(diff(threeLeads(:,3)).^2));
 
end