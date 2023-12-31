%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEMI classifier based on 4th Universal definition of MI j-point
%   criteria
%
% Returns result - Boolean 1 = STEMI detected; 0 = no STEMI detected
% Passed 
%     -averageBeats - N x 9 matrix with one beat. Columns are V1-V6;I-III
%     -jPoint_n - the sample number of the j-point e.g. 460
%     -maleFemale - whether the patient is male (0) or female (1)
%     -age - age in years
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRITERIA: (1 mm == 100uV)
% New ST-elevation at the J-point in two contiguous leads with
% the cut-point: ? 1 mm in all leads other than leads V2–V3 where
% the following cut-points apply: ? 2mm in men ? 40 years;
% ? 2.5 mm in men < 40 years, or ? 1.5 mm in women regardless
% of age
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result, location] = detectSTEMI(averageBeats, jPoint_n, maleFemale, age)

    %expand 9 leads to 12 leads (append aVR, aVL, aVF on the end columns)
%     averageBeats = cell2mat(averageBeats);
    Nleads = length(averageBeats(1,:));
    N = length(averageBeats(:,1));
    beats_12 = zeros(N,12); %preallocation
    beats_12(:,1:9) = averageBeats; 
    beats_12(:,10) = -(0.5)*(beats_12(:,7) + beats_12(:,8));  %aVR = -(1/2)(I + II)
    beats_12(:,11) = beats_12(:,7) - 0.5*beats_12(:,8);  %aVL = I – (1/2) II
    beats_12(:,12) = beats_12(:,8) - beats_12(:,7);  %aVF = II – (1/2)I
    
    %separate the jpoint values into a matrix
    jPoints = zeros(1,12);
    for i = 1:12
        jPoints(i) = beats_12(jPoint_n,i);
    end
   
    %% STEMI CRITERIA %%
    count = 0;
    if(maleFemale == 0) %for males
        %% MORE THAN 40 %%
        if(age >= 40)
            %ANTERIOR (V1-V4)
            for i = 1:4 %for each lead
                if(i==2 || i==3) %V2/V3
                    if(jPoints(i) >= 200e-3) %200uV
                        count = count+1;    %STE present leads +1
                        continue;
                    end
                elseif(jPoints(i) >= 100e-3) %100uV
                    count = count+1;    %STE present leads +1
                    continue;
                end
            end
            if(count > 1)   %if 2 contiguous leads
                result = 1; %STEMI detected
                location = "Anterior";
                return;
            else
                count=0; %else reset the count to 0
            end
            
            %LATERAL (I,aVL,V5,V6)
            lat = [5 6 7 11]; %indexes the above leads in jPoints
            for i = 1:4 %for each lead
                if(jPoints(lat(i)) >= 100e-3) %100uV
                    count = count+1;    %STE present leads +1
                    continue;
                end
            end
            if(count > 1)   %if 2 contiguous leads
                result = 1; %STEMI detected
                location = "Lateral";
                return;
            else
                count=0;    %else reset the count to 0
            end
            
            %INFERIOR (II,III,aVF)
            inf = [8 9 12];
            for i = 1:3 %for each lead
                if(jPoints(inf(i)) >= 100e-3) %100uV
                    count = count+1;    %STE present leads +1
                    continue;
                end
            end
            if(count > 1)   %if 2 contiguous leads
                result = 1; %STEMI detected
                location = "Inferior";
                return;
            else
                count=0;    %else reset the count to 0
            end
  
        %% LESS THAN 40 %%
        else %age < 40 years old
            %ANTERIOR (V1-V4)
            for i = 1:4 %for each lead
                if(i==2 || i==3) %V2/V3
                    if(jPoints(i) >= 250e-3) %200uV
                        count = count+1;
                        continue;
                    end
                elseif(jPoints(i) >= 100e-3) %100uV
                    count = count+1;
                    continue;
                end
            end
            if(count > 1)
                result = 1;
                location = "Anterior";
                return;
            else
                count=0;
            end
            
            %LATERAL (I,aVL,V5,V6)
            lat = [5 6 7 11];
            for i = 1:4 %for each lead
                if(jPoints(lat(i)) >= 100e-3) %100uV
                    count = count+1;
                    continue;
                end
            end
            if(count > 1)
                result = 1;
                location = "Lateral";
                return;
            else
                count=0;
            end
            
            %INFERIOR (II,III,aVF)
            inf = [8 9 12];
            for i = 1:3 %for each lead
                if(jPoints(inf(i)) >= 100e-3) %100uV
                    count = count+1;
                    continue;
                end
            end
            if(count > 1)
                result = 1;
                location = "Inferior";
                return;
            else
                count=0;
            end
        end
    %% FOR FEMALES %%
    else
        %ANTERIOR (V1-V4)
        for i = 1:4 %for each lead
            if(i==2 || i==3) %V2/V3
                if(jPoints(i) >= 150e-3) %200uV
                    count = count+1;
                    continue;
                end
            elseif(jPoints(i) >= 100e-3) %100uV
                count = count+1;
                continue;
            end
        end
        if(count > 1)
            result = 1;
            location = "Anterior";
            return;
        end

        %LATERAL (I,aVL,V5,V6)
        lat = [5 6 7 11];
        for i = 1:4 %for each lead
            if(jPoints(lat(i)) >= 100e-3) %100uV
                count = count+1;
                continue;
            end
        end
        if(count > 1)
            result = 1;
            location = "Lateral";
            return;
        else 
            count=0;
        end

        %INFERIOR (II,III,aVF)
        inf = [8 9 12];
        for i = 1:3 %for each lead
            if(jPoints(inf(i)) >= 100e-3) %100uV
                count = count+1;
                continue;
            end
        end
        if(count > 1)
            result = 1;
            location = "Inferior";
            return;
        else
            count=0;
        end
    end
    %% ST depression
    %ANTERIOR (V1-V4)
    dep = 1000000e-3; %amount of depression required in mV (high number means not required)
    for i = 1:4 %for each lead
        if(jPoints(i) <= -dep) %50uV depression
            count = count+1;
            continue;
        end
    end
    if(count > 1)
        result = 1;
        location = "Anterior Depression";
        return;
    else
        count=0;
    end

    %LATERAL (I,aVL,V5,V6)
    lat = [5 6 7 11];
    for i = 1:4 %for each lead
        if(jPoints(lat(i)) <= -dep) %50uV depression
            count = count+1;
            continue;
        end
    end
    if(count > 1)
        result = 1;
        location = "Lateral Depression";
        return;
    else
        count=0;
    end

    %INFERIOR (II,III,aVF)
    inf = [8 9 12];
    for i = 1:3 %for each lead
        if(jPoints(inf(i)) <= -dep) %50uV depression
            count = count+1;
            continue;
        end
    end
    if(count > 1)
        result = 1;
        location = "Inferior Depression";
        return;
    else
        count=0;
    end
    
    %% Return null
    if(count < 2)   %if no contigous leads detected
        result = 0; %no stemi detected
        location = "null";
        return;
    end
end

%change the J point
%check the ST depression is NOT upsloping

%% PSEUDOCODE
%count=0;
%if male
    %if >40 years old
        %ANTERIOR (V1-V4)
        %for each lead
            %if V2-V3
                %if jpoint > 2mm
                    %count = count+1
            %if jpoint > 1mm
                %count = count+1
        %if count > 1
            %result = 1; location = "Anterior"
            %break;
        %count=0;
        %LATERAL (I,aVL,V5,V6)
        %for each lead
            %if jpoint > 1mm
                %count = count+1
        %if count > 1
            %result = 1; location = "Anterior"
            %break;  
        %count=0;
        %INFERIOR (II,III,aVF)
        %for each lead
            %if jpoint > 1mm
                %count = count+1;
        %if count > 1
            %result = 1; location = "Anterior"
            %break;  
        %count=0
    %if <40yo
        %ANTERIOR (V1-V4)
        %for each lead
            %if V2-V3
                %if jpoint > 2.5mm
                    %count = count+1
            %if jpoint > 1mm
                %count = count+1
        %if count > 1
            %result = 1; location = "Anterior"
            %break;
        %count=0;
        
        %LATERAL (I,aVL,V5,V6)
        %for each lead
            %if jpoint > 1mm
                %count = count+1
        %if count > 1
            %result = 1; location = "Anterior"
            %break;  
        %count=0;
        
        %INFERIOR (II,III,aVF)
        %for each lead
            %if jpoint > 1mm
                %count = count+1;
        %if count > 1
            %result = 1; location = "Anterior"
            %break;  
        %count=0
        
%if female
    %ANTERIOR (V1-V4)
    %for each lead
        %if V2-V3
            %if jpoint > 1.5mm
                %count = count+1
        %if jpoint > 1mm
            %count = count+1
    %if count > 1
        %result = 1; location = "Anterior"
        %break;
    %count=0;
    
    %LATERAL (I,aVL,V5,V6)
    %for each lead
        %if jpoint > 1mm
            %count = count+1
    %if count > 1
        %result = 1; location = "Anterior"
        %break;  
    %count=0;
    
    %INFERIOR (II,III,aVF)
    %for each lead
        %if jpoint > 1mm
            %count = count+1;
    %if count > 1
        %result = 1; location = "Anterior"
        %break;  
    %count=0