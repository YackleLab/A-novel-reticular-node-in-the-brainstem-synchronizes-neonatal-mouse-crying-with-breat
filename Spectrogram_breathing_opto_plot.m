%% This script is to align breathing with spectrogram, and identify/analyze USV
%Pressure Transducer sampling rate of 1000Hz
%Audio sampling rate of 400,000Hz

%% load in breathing data and filtering
pathtodata = '~/Box/Lab/USV Behavior/XYZ';
addpath ('~/Box/Lab/USV Behavior/Code/Whistles');
filename = '_____';
txtfile = strcat (filename,'.txt');
optofile = strcat (filename,'_opto.txt');
cd(pathtodata)
breathdata=readtable(txtfile, 'ReadVariableNames', false);

breathtrace=table2array(breathdata);
breathtracesamprate = 1000.; %hz
time = 0:1/breathtracesamprate:length(breathtrace)/breathtracesamprate;
time = time(1:length(time)-1);

optostim=readtable(optofile, 'ReadVariableNames', false);
optostim=table2array(optostim);
optostim = optostim*1000;
optostim(optostim<2500)=0;
optostim = optostim - 500;

optostim= downsample( optostim ,10 ); 

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

%% Load in audio file and generating spectrogram
[micechirp,fs] = audioread (strcat(filename,'.wav'));
time2 = 0:1/fs:length(micechirp)/fs;
time2 = time2(1:length(time2)-1);
window=flattopwin (512); %window
nooverlap=256; %overlap bw windows
nfft=512; 

[~,F,T,P] = spectrogram (micechirp,window,nooverlap,[],fs,'yaxis');

%% opto sync
% optoall = find (optotrace > 4.9);
% optoend = optoall (diff(optoall)>1);
% optostart = optoend - 10; 
% 
% optoall_extended = find (optotrace_extended > 4.9);
% optoend_extended = optoall_extended (diff(optoall_extended)>1);
% optostart_extended = optoend_extended - 10; 

%% plot spectrogram no tidal volume
figure;
ax1=subplot (2,1,1);
imagesc(T, F, 10*log10(P+eps)); % add eps like pspectrogram does
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (sec)');
% h = colorbar;
% h.Label.String = 'Power/frequency (dB/Hz)';
c = gray;
c = flipud(c);
colormap(c);
caxis([-105 -60]);  

ax2=subplot(2,1,2);
yline(0,':'); 
hold on;
plot (time3, filtered_pressuredif_smooth, 'k');
%   hold on;
%   plot (inspStart/breathtracesamprate, filtered_pressuredif_smooth(inspStart), 'r*');
%   hold on;
%   plot (expStart/breathtracesamprate, filtered_pressuredif_smooth(expStart), 'b*');
hold on;
h = area (time, optostim, -1000, 'LineStyle', 'none');
h.FaceColor = [0, 0.4470, 0.7410];
h.FaceAlpha = 0.1;
ylim([-150 150])
%set(gca,'YTick', []);
xlabel('Time (sec)');
ylabel('Airflow (AU)');
ylim ([-6 6]);

linkaxes ([ax1, ax2], 'x');  

