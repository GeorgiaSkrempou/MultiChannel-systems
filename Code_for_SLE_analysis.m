%% Section 1: read time and channel column
% change filename and range accordingly
time = readtable('2022-01-18T13-47-34BL261_slice1_Recording-0_(Data Acquisition (1);MEA2100-Mini; Electrode Raw Data1)_Analog.csv', 'Range', 'A:A');
voltage = readtable('2022-01-18T13-47-34BL261_slice1_Recording-0_(Data Acquisition (1);MEA2100-Mini; Electrode Raw Data1)_Analog.csv', 'Range', 'R:R');

%% Section 2: Convert units
time_sec=(time{:,1}*1e-6);%converts time from microseconds to seconds
voltage_mV=(voltage{:,1}*1e-9);%converts voltage from picovolt to milivolt

%% Section 3: PLOT LOWPASSED SIGNAL OF THE WHOLE RECORDING
filtered_signal = lowpass(voltage_mV, 150, 10000); %lowpass filter at cutoff 150Hz

figure;
plot(time_sec,voltage_mV, 'color', 'k'); %plots the original signal
xlabel('time (min)');
ylabel('voltage (mV)');
ylim([-1, 1]);
title('Original signal 10kHz');

figure;
plot(time_sec,voltage_mV, time_sec, filtered_signal); %plots the lowpassed signal on top of the original
xlabel('time (min)');
ylabel('voltage (mV)');
ylim([-1, 1]);
title('Superimposed low-passed signal 150Hz');

figure;
plot(time_sec, filtered_signal, 'color', 'k'); %plots the lowpassed signal
xlabel('time (min)');
ylabel('voltage (mV)');
ylim([-1, 1]);
title('low-passed signal 150Hz');

%% Section 4: ANALYZE POWER IN EVERY SLEs
% change the xlsx file directory (leave range unchanged - should be the
% same every time
events = readtable("C:\Users\georg\OneDrive - UvA\Desktop\PhD\MEA project\data\Georgia's analysis\JWH\BL261_events.xlsx", 'Range', 'A:B'); % reads the events from a separate table
signalData = table(time_sec(:), voltage_mV(:), 'VariableNames', {'Time_sec', 'Voltage_mV'}); % makes a table with the time and voltage
rows = height(events); 
Data = array2table(zeros(rows,6));
Data.Properties.VariableNames = {'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma' 'HFOs'};
Fs = 10000; % sampling frequency
N = 3; % filter order
Fn = Fs/2; % Nyiquist frequency
FbpLow = [48 52];
[B, A] = butter(N, [min(FbpLow)/Fn max(FbpLow)/Fn], 'stop'); % butterworth filter


for event = 1:rows

    start_time = table2array(events(event,1));
    end_time = table2array(events(event,2));
    L = table2array(varfun(@(x)((x>=start_time) & (x<=end_time)), signalData(:,1)));
    T2 = signalData(L,:); 
    seizure_time = T2{:,1};
    seizure_voltage = T2{:,2};

    figure;
    plot(seizure_time,seizure_voltage);%Plot SLE
    xlabel('time (min)');
    ylabel('voltage (mV)');
    ylim([-1, 1]);
    title('Seizure-like event');

    %Bandstop filter at [48 52] -filtering the electrical current frequencies
    filt_seizure_voltage_mV = filtfilt(B, A, seizure_voltage);

    % Bandpass filter (Lowpass) at 3000Hz
    filteredLow = lowpass(filt_seizure_voltage_mV, 3000, Fs);
    %figure()
    %plot(seizure_time,filt_seizure_voltage_mV)

    % Power-Spectrum analysis and plotting
    signalPS = filteredLow; 
    L = length(signalPS);
    fftSignal = fft(signalPS);
    p2 = abs(fftSignal/L);
    p1 = p2(1:L/2+1);
    p1(2:end-1) = 2*p1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    figure()
    plot(f,p1) %plot power
    title('Single-Sided Amplitude Spectrum of signal')
    xlabel('f (Hz)')
    ylabel('|p1(f)|')


    % CALCULATE POWER 
    Delta = bandpower(filteredLow,Fs,[1 4]);
    Theta = bandpower(filteredLow,Fs,[4 8]);
    Alpha = bandpower(filteredLow,Fs,[8 12]);
    Beta = bandpower(filteredLow,Fs,[12 30]);
    Gamma = bandpower(filteredLow,Fs,[30 100]);
    HFOs = bandpower(filteredLow,Fs,[100 600]);

    tempData = [Delta,Theta,Alpha,Beta,Gamma, HFOs];

    Data(event,1) = {Delta};
    Data(event,2) = {Theta};
    Data(event,3) = {Alpha};
    Data(event,4) = {Beta};
    Data(event,5) = {Gamma};
    Data(event,6) = {HFOs};
  
    figure;
    name = {'Delta';'Theta';'Alpha';'Beta';'Gamma';'HFOs'};
    x = [1:6]; y = tempData; 
    bar(x,y)
    set(gca,'xticklabel',name)
    xlabel('Frequency (Hz)');
    ylabel('Power (mV^2)');
    title('Event powers of frequency bands')

end
%% Section 5: HIGHPASS FILTER TO SEE NEURONAL FIRING 
%filteredHigh = bandpass(seizure_voltage,[300 3000],Fs);
%figure;
%plot(seizure_time,filteredHigh);
%xlabel('time (min)');
%ylabel('voltage (mV)');
%ylim([-1, 1]);
%title('Seizure-like event');


