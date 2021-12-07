%% This script is to align breathing with spectrogram, and identify/analyze USV
%Pressure Transducer sampling rate of 1000Hz
%Audio sampling rate of 400,000Hz

%% load in breathing data and filtering
tic
condition = 'XYZ';
pathtodata = strcat('~/Box/Lab/USV Behavior/', condition);
condition = strrep(condition, '/', ' ');
condition = erase (condition, ' grouped');


addpath ('~/Box/Lab/USV Behavior/Code/Whistles')

aggbasalbreathDur = [];
aggboutbreathDur = [];

d=dir(fullfile(pathtodata,'*.wav'));
file_names={d.name};
for k=1:numel(file_names)
 try 

%load in breathing data
[~, filename, ~] = fileparts(file_names{k});

%if isfile (strcat(filename,'breathing_vars_072121.mat')) == 0
txtfile = strcat (filename,'.txt');
cd(pathtodata)
breathdata=readtable(txtfile, 'ReadVariableNames', false);
breathtrace=table2array(breathdata);
breathtracesamprate = 1000.; %hz
time = 0:1/breathtracesamprate:length(breathtrace)/breathtracesamprate;
time = time(1:length(time)-1);

%Bandpass filtering data at 2Hz to 35hz - can adjust frequency if baseline looks
%really off
filtered_breathtrace=bandpass(breathtrace,[2,35],1000);%

%airflow using filtered tidal volume
pressuredif = minus (filtered_breathtrace (1:length(filtered_breathtrace)-1), filtered_breathtrace(2:length(filtered_breathtrace))); 
pressuredif = pressuredif*1000;
time3  = time (1:length(time)-1);
filtered_pressuredif = highpass(pressuredif,2,1000); 
filtered_pressuredif_smooth = smoothdata (filtered_pressuredif,'movmean',30);

%% Getting Breath Parameters
%find min and max of tidal volume (aka start and end of breath)
[maxpks, ~] = findpeaks(filtered_breathtrace ,'MinPeakProminence',0.030);
[minpks, ~] = findpeaks(-1*filtered_breathtrace ,'MinPeakProminence',0.030);
MPP = (median(minpks)+median(maxpks))/2.5;
[maxpks, localmax] = findpeaks(filtered_breathtrace ,'MinPeakProminence',MPP);
[minpks, localmin] = findpeaks(-1*filtered_breathtrace ,'MinPeakProminence',MPP);

%make sure it starts with onset of inspiration and ends with offset of expiration   
if localmax (1) < localmin (1)
    localmax = localmax (2:end);
end
if localmax(end)>localmin(end) 
    localmax = localmax (1:end-1);
end

%define inspiratory start time, expiratory start time and expiratory end
%time (should all have same number of elements)
inspStart = localmin(1:end-1);
expStart = localmax;
expEnd = localmin(2:end);
inspDur = expStart-inspStart;
expDur = expEnd-expStart;
breathDur = expEnd - inspStart; 

%thresholding on inspiratory duration and breath duration
%minBreathDurThresh = 110;
inspStart = inspStart(inspDur > 40);
expStart = expStart(inspDur > 40);
expEnd = vertcat(inspStart(2:end), expEnd(end));
%recalculate breathparameters after  
inspDur = expStart-inspStart;
expDur = expEnd - expStart; 
breathDur = expEnd - inspStart; 
instfreq = 1000./breathDur;
tidvolpks = filtered_breathtrace (expStart);

% Filing inspiratory and expiratory peak values
inspPeak = [];
expPeak = [];
for i= 1 : length (inspStart)
    breath = filtered_pressuredif_smooth (inspStart(i):expEnd(i)-1);
    inspPeak = [inspPeak min(breath)];
    expPeak = [expPeak max(breath)];
end

   [micechirp,fs] = audioread (strcat(filename,'.wav'));
    time2 = 0:1/fs:length(micechirp)/fs;
    time2 = time2(1:length(time2)-1);
   
    
   %% Run Holy Lab Tools to identify USVs
    %initialize Parameters for Sound2Sng
    sngparms.plot = false;
    sngparms.threshold = 1010;%900;
    sngparms.nfreq = 256;
    lowbound=10000;
    upperbound=150000;
    sngparms.freqrange = [lowbound upperbound];

    %initialize Parameters for Whistimes
    whistimesparms.puritythresh = 0.3;
    whistimesparms.specdiscthresh = 0.8;
    whistimesparms.durationthresh = 0.002;
    whistimesparms.mergeclose = 0.015;
    whistimesparms.meanfreqthresh = 30000;

    %songname
    sngname = strcat ('sng_',filename);

    if isfile (sngname) == 0
        %generate sngname and micechirp (voltage trace)
        sound2sng(strcat (filename,'.wav'),sngparms,sngname);
    else
    end 

    %get chirp event times
    twhis = whistimes(sngname,whistimesparms);

    twhis2 = twhis*fs;




    %% Timing of Vocalization relative to breathing parameters

    % %find values closest and previous to set time (aka last inspiration that
    % %occured before a vocalization)
    vocStartTime = twhis (1,:)*1000;
    vocEndTime = twhis (2,:)*1000;
    if vocStartTime(1)<inspStart(1)
    vocStartTime = vocStartTime(2:end);
    vocEndTime = vocEndTime(2:end);
    end 
       
    vocBreathInspStartInd = interp1 (inspStart, 1:numel(inspStart), vocStartTime, 'previous', 'extrap');
    vocExpStart = transpose(expStart(vocBreathInspStartInd));

    vocTVpks = sort (filtered_breathtrace(vocExpStart));
    minVocTV = max (0.05, prctile (vocTVpks, 5));




      
  breathDur = breathDur.';
  basalbreathDur = breathDur (tidvolpks<0.05);
  basalbreathDur2 = prctile (basalbreathDur, [5,95]);
  basalbreathDur3 = basalbreathDur(basalbreathDur>basalbreathDur2(1) & basalbreathDur<basalbreathDur2(2));
  
  boutbreathDur = breathDur (tidvolpks>minVocTV);
  boutbreathDur2 = prctile (boutbreathDur, [5,95]);
  boutbreathDur3 = boutbreathDur(boutbreathDur>boutbreathDur2(1) & boutbreathDur<boutbreathDur2(2));
  
%   basalbreathfreq = 1./((breathDur (tidvolpks<0.05))/1000.); 
%   basalbreathfreq2 = prctile (basalbreathfreq, [5,95]);
%   basalbreathfreq3 = basalbreathfreq(basalbreathfreq>basalbreathfreq2(1) & basalbreathfreq<basalbreathfreq2(2));
  
 aggbasalbreathDur = [aggbasalbreathDur basalbreathDur3];
 aggboutbreathDur = [aggboutbreathDur boutbreathDur3];
 
 
 
  catch
       fprintf('loop number %d failed\n',k)
end
 
end  
%%
%  figure; 
%  xhistfit (aggbreathfreq);
%  figure;
%  histfit (aggbreathDur);

aggbasalbreathDur = aggbasalbreathDur(aggbasalbreathDur<1000);
 
aggboutbreathDur = aggboutbreathDur (aggboutbreathDur <1000);

basalmatrix = [(aggbasalbreathDur(1:length(aggbasalbreathDur)-1));(aggbasalbreathDur(2:length(aggbasalbreathDur)))];

meanbasaln = mean(basalmatrix(1,:));
%meanbasaln1 = mean(basalmatrix(2,:));
stdbasaln = std(basalmatrix(1,:));
%stdbasaln1 = std(basalmatrix(2,:));

% Create rotation matrix
theta = -45; % to rotate 90 counterclockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% Rotate your point(s)
basalrotpoint = R*basalmatrix;

stdbasallong = std(basalrotpoint(1,:));
stdbasalshort = std(basalrotpoint(2,:));

boutmatrix = [(aggboutbreathDur(1:length(aggboutbreathDur)-1));(aggboutbreathDur(2:length(aggboutbreathDur)))];

meanboutn = mean(boutmatrix(1,:));
stdboutn = std(boutmatrix(1,:));

% Create rotation matrix
theta = -45; % to rotate 90 counterclockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% Rotate your point(s)
boutrotpoint = R*boutmatrix;

stdboutlong = std(boutrotpoint(1,:));
stdboutshort = std(boutrotpoint(2,:));



figure;
 plot (aggbasalbreathDur (1:length(aggbasalbreathDur)-1), aggbasalbreathDur (2:length(aggbasalbreathDur)),'.');
 xlabel({'Breath duration (n)','(ms)'})
 ylabel({'Breath duration (n+1)','(ms)'})
 title(strcat(condition, ' Basal Breath Durations'))
 xlim([0 1000])
 ylim([0 1000]) 
 hold on;
 refline (1,0);
pbaspect([1 1 1])
hold on;
errorbar(meanbasaln,meanbasaln,stdbasaln,'both')

basallongx = [meanbasaln-(stdbasallong)/sqrt(2) meanbasaln+(stdbasallong)/sqrt(2)];
basallongy = [meanbasaln-(stdbasallong)/sqrt(2) meanbasaln+(stdbasallong)/sqrt(2)];
basalshortx = [meanbasaln-(stdbasalshort)/sqrt(2) meanbasaln+(stdbasalshort)/sqrt(2)];
basalshorty = [meanbasaln+(stdbasalshort)/sqrt(2) meanbasaln-(stdbasalshort)/sqrt(2)];

hold on;
line(basallongx,basallongy,'Color','#D95319','Marker','x')
hold on;
line(basalshortx,basalshorty,'Color','#D95319','Marker','x')
saveas(gcf,strcat(condition,' Basal Breath Durations'), 'epsc');
 

figure;
 plot (aggboutbreathDur (1:length(aggboutbreathDur)-1), aggboutbreathDur (2:length(aggboutbreathDur)),'.');
 xlabel({'Breath duration (n)','(ms)'})
 ylabel({'Breath duration (n+1)','(ms)'})
 title(strcat(condition,' Bout Breath Durations'))
 xlim([0 1000])
 ylim([0 1000]) 
 hold on;
 refline (1,0);
pbaspect([1 1 1])

errorbar(meanboutn,meanboutn,stdboutn,'both')

boutlongx = [meanboutn-(stdboutlong)/sqrt(2) meanboutn+(stdboutlong)/sqrt(2)];
boutlongy = [meanboutn-(stdboutlong)/sqrt(2) meanboutn+(stdboutlong)/sqrt(2)];
boutshortx = [meanboutn-(stdboutshort)/sqrt(2) meanboutn+(stdboutshort)/sqrt(2)];
boutshorty = [meanboutn+(stdboutshort)/sqrt(2) meanboutn-(stdboutshort)/sqrt(2)];

hold on;
line(boutlongx,boutlongy,'Color','#D95319','Marker','x')
hold on;
line(boutshortx,boutshorty,'Color','#D95319','Marker','x')
saveas(gcf,strcat (condition,' Bout Breath Durations'), 'epsc');




toc

%end 

  
    