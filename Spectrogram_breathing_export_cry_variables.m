%% This script is to align breathing with spectrogram, and identify/analyze USV
%Pressure Transducer sampling rate of 1000Hz
%Audio sampling rate of 400,000Hz

%% load in breathing data and filtering
pathtodata = '~/Box/Lab/USV Behavior/_____';
addpath ( '~/Box/Lab/USV Behavior/Code/Whistles')

d=dir(fullfile(pathtodata,'*.wav'));
file_names={d.name};
for k=1:numel(file_names)
 try 
tic
%load in breathing data
[~, filename, ~] = fileparts(file_names{k});

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
%alternative to use 2nd order butterworth filter
%[butterb, buttera] = butter (2, [2 35]/500);
%filtered_breathtrace2=filtfilt(butterb,buttera,breathtrace); 

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
%%%%% plotting some initial distribution of breathing parameter 


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


%% Load in audio file and generating spectrogram
[micechirp,fs] = audioread (strcat(filename,'.wav'));
time2 = 0:1/fs:length(micechirp)/fs;
time2 = time2(1:length(time2)-1);
window=flattopwin (512); %window
nooverlap=256; %overlap bw windows
nfft=512; 

[~,F,T,P] = spectrogram (micechirp,window,nooverlap,[],fs,'yaxis');


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
%nextVocBreathInspStartInd = interp1 (inspStart, 1:numel(inspStart), vocEndTime, 'next', 'extrap');
%nextVocInspStart = transpose(inspStart(nextVocBreathInspStartInd));

vocInspStart = transpose(inspStart(vocBreathInspStartInd));
vocBreathDur = transpose (breathDur(vocBreathInspStartInd));
vocExpDur = transpose (expDur(vocBreathInspStartInd));
delayOnsetExp = vocStartTime - vocExpStart;
delayOffsetExp = vocEndTime-vocExpStart;
delayOnsetInsp = vocStartTime - vocInspStart;
delayOffsetInsp = vocEndTime-vocInspStart;
normDelayOnsetExp = delayOnsetExp./vocExpDur;
normDelayOffsetExp =  delayOffsetExp./vocExpDur;
normDelayOnsetInsp = delayOnsetInsp./vocBreathDur;
normDelayOffsetInsp = delayOffsetExp./vocBreathDur;

multisyllabindex = or ([diff(vocExpStart) 1] == 0, [1 diff(vocExpStart)] == 0);
singlesyllabindex = ~multisyllabindex;
[uniqcalls, ~, rnk] = unique(vocInspStart);
rnk = transpose (rnk);

numcalls = length (uniqcalls); 
numMScalls = sum(histc(vocInspStart, unique(vocInspStart)) > 1);

[callscounts, ~] = histc(vocInspStart, uniqcalls);
bisyllabvocstarts = uniqcalls(callscounts ==2);
trisyllabvocstarts = uniqcalls(callscounts >=3);
[bisyllabindex,~] = ismember(vocInspStart,bisyllabvocstarts);
[trisyllabindex, ~] = ismember(vocInspStart,trisyllabvocstarts);

unisyllabicdelayOnsetExp = delayOnsetExp(singlesyllabindex);
unisyllabicNormDelayOnsetExp = normDelayOnsetExp(singlesyllabindex);
multiDelayOnsetExp = delayOnsetExp(multisyllabindex);
multiNormDelayOnsetExp = normDelayOnsetExp(multisyllabindex);
biDelayOnsetExp = delayOnsetExp(bisyllabindex);
biNormDelayOnsetExp = normDelayOnsetExp(bisyllabindex);
triDelayOnsetExp = delayOnsetExp(trisyllabindex);
triNormDelayOnsetExp = normDelayOnsetExp(trisyllabindex);

save (fullfile(pathtodata,strcat(filename,'usv_culmative_vars.mat')),'numcalls','numMScalls', 'delayOnsetExp','normDelayOnsetExp', 'unisyllabicdelayOnsetExp', 'unisyllabicNormDelayOnsetExp','multiDelayOnsetExp','multiNormDelayOnsetExp','biDelayOnsetExp','biNormDelayOnsetExp','triDelayOnsetExp','triNormDelayOnsetExp');

toc
 catch
       fprintf('loop number %d failed\n',k)
 end
end 
