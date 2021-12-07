clear
tic
filename = strcat(date,'_opto_power_hz_sec_');
s = daq.createSession('ni');  
s.Rate = 250000;
s.DurationInSeconds = 120;

channel1 = addAnalogInputChannel(s,'Dev1','ai0','Voltage');
channel1.Range = [-5.000000 5.000000];
 
channel2 = addAnalogInputChannel(s,'Dev1','ai1','Voltage');
channel2.TerminalConfig = 'SingleEnded';
channel2.Range = [-5.0000000 5.0000000];

channel3 = addAnalogInputChannel(s,'Dev1','ai2','Voltage');
channel3.TerminalConfig = 'SingleEnded';
channel3.Range = [-10.0000000 10.0000000];

[data, timestamps, starttime] = startForeground(s);
audiodata = data (:,1);
airflowdata = data (:,2);
optodata = data (:,3);
airflowdata = downsample (airflowdata, s.Rate/1000);
writematrix(airflowdata, strcat ('C:\Users\Yackle Lab\Documents\USV Data\',filename,'.txt'));

optodata = downsample (optodata, s.Rate/1000);
writematrix(optodata, strcat ('C:\Users\Yackle Lab\Documents\USV Data\',filename,'_opto.txt'));


audiodata_normalized = audiodata/5; 
audiowrite (strcat ('C:\Users\Yackle Lab\Documents\USV Data\',filename,'.wav'), audiodata_normalized, s.Rate);
beep;
 