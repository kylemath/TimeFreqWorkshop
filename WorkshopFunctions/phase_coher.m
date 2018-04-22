function [all_coher,all_phasediff] = phase_coher(Y,Fs,window_size,hipass,lopass)
%% phase_coher.m - Created Feb 14, 2012 by Kyle Mathewson, University of Illinois 
% Takes the time x channel data matrix Y, filters out a single frequency
% band, and computes the phase coherence between each channel
% This is done over time by considering the data in discrete windows
% Output is all_coher (the coherence matrix) and all_phasediff (the average phase differences)
% Input is: Y  - data matrix, time points by channels
%           Fs - Sampling Frequency in Hz
%           window_size - size of analysis window in seconds
%           hipass - lower bound of filter
%           lopass - upper bound of filter


%--------------------------------------------------------------------------
%% Initialize some variables
Y_slow = Y;    %Y would be your data file with time points by channels
n_chan = size(Y,2);  %number of channels
step = Fs*window_size; %Size of each coherence window to analyze in seconds;
n_steps = floor(length(Y)/step); %Count how many windows overall
all_coher = zeros(n_steps,n_chan,n_chan);  %Allocate memory for output
all_phasediff = all_coher;

fprintf(['The data is ' num2str(length(Y)/Fs) ' seconds long. ' '\n']); 
fprintf(['Window size is ' num2str(window_size) ' seconds long. ' '\n']);



%--------------------------------------------------------------------------
%% Filter to get the Alpha
order = 3;
bandpass = [(hipass*2)/Fs, (lopass*2)/Fs];          % Bandwidth of bandpass filter
[Bbp,Abp] = butter(order,bandpass);                 % Generation of Xth order Butterworth highpass filter
for c = 2:size(Y,2)
    Y_slow(:,c) = filtfilt(Bbp,Abp,(Y(:,c)));       % Butterworth bandpass filtering of YY
end



%-------------------------------------------------------------------------
%% Compute the alpha coherence between channels and compute the MDS solution
Y_slow_phase = angle(hilbert(Y_slow));              % Take the instantaneous angle of the filtered time series in radians using a hilbert transform

frame_count = 0;  %Counter for each time step
for i_time = 1:step:length(Y)-step                  % for each time window
    tic                                             % start a clock
    range = i_time:i_time+step;                     % find the time range based on the current window
    frame_count = frame_count + 1;                  % increment counter
    fprintf(['Processing window ' num2str(frame_count) ' of ' num2str(n_steps) ' windows.  ']);
    
    for i_P1 = 1:n_chan                             % cycle through each channel and take out the data
        seg_P1 = Y_slow_phase(range,i_P1);          % data segment
    
        for i_P2 = 1:n_chan                         % cycle for the channel with which to compare  
            if i_P2 == i_P1                         % in order to avoid rounding errors 
                all_coher(frame_count,i_P1,i_P2) = 1;
                all_phasediff(frame_count,i_P1,i_P2) = 0;
            else                                    % Compute the phase coherence and record the values
                seg_P2 = Y_slow_phase(range,i_P2);  % data segment
                mean_vec = sum(exp(1i*(seg_P1-seg_P2)))/step;  %compute the phase diff  in imaginary plane as the length of the average phase difference vector 
                all_coher(frame_count,i_P1,i_P2) = abs(mean_vec);  %record the values
                all_phasediff(frame_count,i_P1,i_P2) = angle(mean_vec);
            end
        end
        
    end
    
    toc                                             %stop clock and display elapsed time
end
