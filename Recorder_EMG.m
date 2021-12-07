clear
filename = strcat(date,'emg_cage_x_mouse_y');
s = daq.createSession('ni');  
s.Rate = 250000;
s.DurationInSeconds = 120;

channel1 = addAnalogInputChannel(s,'Dev1','ai0','Voltage');
channel1.Range = [-5.000000 5.000000];
 
channel2 = addAnalogInputChannel(s,'Dev1','ai1','Voltage');
channel2.TerminalConfig = 'SingleEnded';
channel2.Range = [-1.0000000 1.0000000];

channel3 = addAnalogInputChannel(s,'Dev1','ai2','Voltage');
channel3.TerminalConfig = 'SingleEnded';
channel3.Range = [-5.0000000 5.0000000];

[data, timestamps, starttime] = startForeground(s);
audiodata = data (:,1);
airflowdata = data (:,2);
emgdata = data (:,3);
airflowdata = downsample (airflowdata, s.Rate/1000);
writematrix(airflowdata, strcat ('C:\Users\Yackle Lab\Documents\USV Data\',filename,'.txt'));

emgdata = downsample (emgdata, s.Rate/10000);
writematrix(emgdata, strcat ('C:\Users\Yackle Lab\Documents\USV Data\',filename,'_emg.txt'));


audiodata_normalized = audiodata/5; 
audiowrite (strcat ('C:\Users\Yackle Lab\Documents\USV Data\',filename,'.wav'), audiodata_normalized, s.Rate);
beep;
  
