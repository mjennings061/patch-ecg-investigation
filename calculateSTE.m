%% Calculate the mean and median ST-elevation across the 12-lead ECG, Horacek VSLs and a chosen short spaced lead
%Input:
% - BalloonBSPMdata - 352 node BSPM data
% - patientNo - Chosen patient index numbers for testing
% - jDelay - the amount of samples after the J point for STE calculation
%Output:
% - avg_median - a matrix of mean and median STE for each lead

function avg_median = calculateSTE(BalloonBSPMdata, patientNo, jDelay, node_p, node_n)
tic
    avg_median = [];
    %% Loop through all positively responding patients
    for h = 1:2:length(patientNo)-1
        %% Find j-point and amplitude for baseline and peak leads
        for i=1:length(BalloonBSPMdata{patientNo(h)})    %for all baseline annotations
            if(BalloonBSPMdata{patientNo(h)}(1,i) == 3)        %find the first j point annotation
                jPoint = i+jDelay;          %add the 40ms delay
                amp_baseline_svl = BalloonBSPMdata{patientNo(h)}(node_p+3,jPoint) - BalloonBSPMdata{patientNo(h)}(node_n+3,jPoint);
                amp_baseline_RA = (BalloonBSPMdata{patientNo(h)}(60+3,jPoint) + BalloonBSPMdata{patientNo(h)}(101+3,jPoint))/2;
                amp_baseline_LA = (BalloonBSPMdata{patientNo(h)}(50+3,jPoint) + BalloonBSPMdata{patientNo(h)}(90+3,jPoint))/2;
                amp_baseline_LL = (3*BalloonBSPMdata{patientNo(h)}(343+3,jPoint) + 2*BalloonBSPMdata{patientNo(h)}(344+3,jPoint))/5;
                amp_baseline_I = amp_baseline_LA - amp_baseline_RA;
                amp_baseline_II = amp_baseline_LL - amp_baseline_RA;
                amp_baseline_III = amp_baseline_LL - amp_baseline_LA;
                amp_baseline_aVR = -(amp_baseline_I + amp_baseline_II)/2;
                amp_baseline_aVL = (amp_baseline_I - amp_baseline_III)/2;
                amp_baseline_aVF = (amp_baseline_II + amp_baseline_III)/2;
                amp_baseline_v1 = BalloonBSPMdata{patientNo(h)}(169+3,jPoint);
                amp_baseline_v2 = BalloonBSPMdata{patientNo(h)}(171+3,jPoint);
                amp_baseline_v3 = (BalloonBSPMdata{patientNo(h)}(192+3,jPoint) + BalloonBSPMdata{patientNo(h)}(193+3,jPoint))/2;
                amp_baseline_v4 = BalloonBSPMdata{patientNo(h)}(216+3,jPoint);
                amp_baseline_v5 = (BalloonBSPMdata{patientNo(h)}(217+3,jPoint) + 2*BalloonBSPMdata{patientNo(h)}(218+3,jPoint))/3;
                amp_baseline_v6 = BalloonBSPMdata{patientNo(h)}(219+3,jPoint);
                amp_baseline_vsl_LAD = BalloonBSPMdata{patientNo(h)}(174+3,jPoint) - BalloonBSPMdata{patientNo(h)}(221+3,jPoint);
                amp_baseline_vsl_LCX = BalloonBSPMdata{patientNo(h)}(221+3,jPoint) - BalloonBSPMdata{patientNo(h)}(150+3,jPoint);
                amp_baseline_vsl_RCA = BalloonBSPMdata{patientNo(h)}(342+3,jPoint) - BalloonBSPMdata{patientNo(h)}(129+3,jPoint);
                break
            end
        end
        for i=1:length(BalloonBSPMdata{patientNo(h)+1})    %for all baseline annotations
            if(BalloonBSPMdata{patientNo(h)+1}(1,i) == 3)        %find the first j point annotation
                jPoint_2 = i+jDelay;          %add the 40ms delay
                amp_peak_svl = BalloonBSPMdata{patientNo(h)+1}(node_p+3,jPoint_2)-BalloonBSPMdata{patientNo(h)+1}(node_n+3,jPoint_2);
                amp_peak_RA = (BalloonBSPMdata{patientNo(h)+1}(60+3,jPoint_2) + BalloonBSPMdata{patientNo(h)+1}(101+3,jPoint_2))/2;
                amp_peak_LA = (BalloonBSPMdata{patientNo(h)+1}(50+3,jPoint_2) + BalloonBSPMdata{patientNo(h)+1}(90+3,jPoint_2))/2;
                amp_peak_LL = (3*BalloonBSPMdata{patientNo(h)+1}(343+3,jPoint_2) + 2*BalloonBSPMdata{patientNo(h)+1}(344+3,jPoint_2))/5;
                amp_peak_I = amp_peak_LA - amp_peak_RA;
                amp_peak_II = amp_peak_LL - amp_peak_RA;
                amp_peak_III = amp_peak_LL - amp_peak_LA;
                amp_peak_aVR = -(amp_peak_I + amp_peak_II)/2;
                amp_peak_aVL = (amp_peak_I - amp_peak_III)/2;
                amp_peak_aVF = (amp_peak_II + amp_peak_III)/2;
                amp_peak_v1 = BalloonBSPMdata{patientNo(h)+1}(169+3,jPoint_2);
                amp_peak_v2 = BalloonBSPMdata{patientNo(h)+1}(171+3,jPoint_2);
                amp_peak_v3 = (BalloonBSPMdata{patientNo(h)+1}(192+3,jPoint_2) + BalloonBSPMdata{patientNo(h)+1}(193+3,jPoint_2))/2;
                amp_peak_v4 = BalloonBSPMdata{patientNo(h)+1}(216+3,jPoint_2);
                amp_peak_v5 = (BalloonBSPMdata{patientNo(h)+1}(217+3,jPoint_2) + 2*BalloonBSPMdata{patientNo(h)+1}(218+3,jPoint_2))/3;
                amp_peak_v6 = BalloonBSPMdata{patientNo(h)+1}(219+3,jPoint_2);
                amp_peak_vsl_LAD = BalloonBSPMdata{patientNo(h)+1}(174+3,jPoint_2) - BalloonBSPMdata{patientNo(h)+1}(221+3,jPoint_2);
                amp_peak_vsl_LCX = BalloonBSPMdata{patientNo(h)+1}(221+3,jPoint_2) - BalloonBSPMdata{patientNo(h)+1}(150+3,jPoint_2);
                amp_peak_vsl_RCA = BalloonBSPMdata{patientNo(h)+1}(342+3,jPoint_2) - BalloonBSPMdata{patientNo(h)+1}(129+3,jPoint_2);
                break
            end
        end
        delta_svl = amp_peak_svl - amp_baseline_svl; %ST elevation from base to peak
        delta_I = amp_peak_I - amp_baseline_I;
        delta_II = amp_peak_II - amp_baseline_II;
        delta_III = amp_peak_III - amp_baseline_III;
        delta_aVR = amp_peak_aVR - amp_baseline_aVR;
        delta_aVL = amp_peak_aVL - amp_baseline_aVL;
        delta_aVF = amp_peak_aVF - amp_baseline_aVF;
        delta_v1 = amp_peak_v1 - amp_baseline_v1;
        delta_v2 = amp_peak_v2 - amp_baseline_v2;
        delta_v3 = amp_peak_v3 - amp_baseline_v3;
        delta_v4 = amp_peak_v4 - amp_baseline_v4;
        delta_v5 = amp_peak_v5 - amp_baseline_v5;
        delta_v6 = amp_peak_v6 - amp_baseline_v6;
        delta_vsl_LAD = amp_peak_vsl_LAD - amp_baseline_vsl_LAD; %three Horacek VSLs
        delta_vsl_LCX = amp_peak_vsl_LCX - amp_baseline_vsl_LCX;
        delta_vsl_RCA = amp_peak_vsl_RCA - amp_baseline_vsl_RCA;

        %store all calculations in a new row of an array
        avg_median = [ avg_median; patientNo(h) delta_svl delta_I delta_II ...
                        delta_III delta_aVR delta_aVL delta_aVF delta_v1 ...
                        delta_v2 delta_v3 delta_v4 delta_v5 delta_v6 ...
                        delta_vsl_LAD delta_vsl_LCX delta_vsl_RCA]; %add the amp_delta_lead... here
    end

    %work out the mean, median outside the loop here
    avg_svl = mean(abs(avg_median(:,2)));
    avg_I = mean(abs(avg_median(:,3)));
    avg_II = mean(abs(avg_median(:,4)));
    avg_III = mean(abs(avg_median(:,5)));
    avg_aVR = mean(abs(avg_median(:,6)));
    avg_aVL = mean(abs(avg_median(:,7)));
    avg_aVF = mean(abs(avg_median(:,8)));
    avg_v1 = mean(abs(avg_median(:,9)));
    avg_v2 = mean(abs(avg_median(:,10)));
    avg_v3 = mean(abs(avg_median(:,11)));
    avg_v4 = mean(abs(avg_median(:,12)));
    avg_v5 = mean(abs(avg_median(:,13)));
    avg_v6 = mean(abs(avg_median(:,14)));
    avg_vsl_LAD = mean(abs(avg_median(:,15)));
    avg_vsl_LCX= mean(abs(avg_median(:,16)));
    avg_vsl_RCA = mean(abs(avg_median(:,17)));

    median_svl = median(abs(avg_median(:,2)));
    median_I = median(abs(avg_median(:,3)));
    median_II = median(abs(avg_median(:,4)));
    median_III = median(abs(avg_median(:,5)));
    median_aVR = median(abs(avg_median(:,6)));
    median_aVL = median(abs(avg_median(:,7)));
    median_aVF = median(abs(avg_median(:,8)));
    median_v1 = median(abs(avg_median(:,9)));
    median_v2 = median(abs(avg_median(:,10)));
    median_v3 = median(abs(avg_median(:,11)));
    median_v4 = median(abs(avg_median(:,12)));
    median_v5 = median(abs(avg_median(:,13)));
    median_v6 = median(abs(avg_median(:,14)));
    median_vsl_LAD = median(abs(avg_median(:,15)));
    median_vsl_LCX = median(abs(avg_median(:,16)));
    median_vsl_RCA = median(abs(avg_median(:,17)));

    %% Put all averages and medians in two new rows and clear unused vars
    avg_median = [avg_median; 1111 avg_svl avg_I avg_II avg_III avg_aVR ...
                                avg_aVL avg_aVF avg_v1 avg_v2 avg_v3 ...
                                avg_v4 avg_v5 avg_v6 avg_vsl_LAD ...
                                avg_vsl_LCX avg_vsl_RCA; ...
                                5050 median_svl median_I median_II median_III ...
                                median_aVR median_aVL median_aVF ...
                                median_v1 median_v2 median_v3 median_v4 ...
                                median_v5 median_v6 median_vsl_LAD ...
                                median_vsl_LCX median_vsl_RCA];
t = toc;
disp(['calculateSTE: ', num2str(t), ' seconds']);
end