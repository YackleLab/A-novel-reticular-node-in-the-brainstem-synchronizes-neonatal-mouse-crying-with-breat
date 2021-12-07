%% This script is to align breathing with spectrogram, and identify/analyze USV
%Pressure Transducer sampling rate of 1000Hz

%% load in breathing data and filtering
pathtodata = '~/Box/Lab/USV Behavior/';
addpath ('~/Box/Lab/USV Behavior/Code/Whistles')
filename = 'abc';
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
%figure;
% nhist(inspDur,'samebins','minbins',25,'smooth');
% xlim ([0 300]); 
% mean (inspDur)

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

%% Plot spectrogram
%plot spectrogram
figure;
ax1=subplot (3,1,1);
imagesc(T, F, 10*log10(P+eps)); % add eps like pspectrogram does
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (sec)');
% h = colorbar;
% h.Label.String = 'Power/frequency (dB/Hz)';
caxis([-105 -60]);


ax2=subplot(3,1,2);

plot(time,filtered_breathtrace,'k');
hold on;
yline(0,':'); 
% hold on;
%  plot (inspStart/breathtracesamprate, filtered_breathtrace(inspStart), 'r*');
%  hold on;
%  plot (expStart/breathtracesamprate, filtered_breathtrace(expStart), 'b*');
 hold on;
  plot([twhis(1,:);twhis(2,:)], [0;0],'LineWidth',5, 'Color', [0.4660 0.6740 0.1880]);
  hold on;
set(gca,'YTick', []);
xlabel('Time (sec)');
ylabel('Tidal Volume');


ax3=subplot(3,1,3);
% plot(time3,filtered_pressuredif,'k');
% hold on;

%%multiplet identification

%[maxintpeaks, intpeakstime] = findpeaks(filtered_pressuredif_smooth ,'MinPeakProminence',0.1);
%[minintpeaks, minintpeakstime] = findpeaks(filtered_pressuredif_smooth * -1 ,'MinPeakProminence',4);


yline(0,':'); 
hold on;
plot (time3, filtered_pressuredif_smooth, 'k');
%  hold on;
%  plot (inspStart/breathtracesamprate, filtered_pressuredif_smooth(inspStart), 'r*');
%  hold on;
%  plot (expStart/breathtracesamprate, filtered_pressuredif_smooth(expStart), 'b*');
%hold on;
% plot (intpeakstime/breathtracesamprate, maxintpeaks, 'g*');
%  hold on;
% plot (minintpeakstime/breathtracesamprate, minintpeaks*-1, 'k*');
 
  hold on;
  plot([twhis(1,:);twhis(2,:)], [0,0],'LineWidth',5, 'Color', [0.4660 0.6740 0.1880]);
ylabel('Airflow') ;
xlabel('Time (sec)');

% Plot 
% ax4=subplot(4,1,4);
% plot (inspStart/breathtracesamprate, instfreq, 'g*');
% xlabel('Time (sec)') ;
% ylabel('Instaneous Freq (Hz)');

linkaxes ([ax1, ax2, ax3], 'x');  

