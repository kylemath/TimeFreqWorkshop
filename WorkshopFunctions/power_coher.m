function [all_coher] = power_coher(Y,Fs,window_size,hipass,lopass,numboot)
%% power_coher.m - Created Feb 21, 2012 by Kyle Mathewson, University of Illinois 
%[all_coher] = power_coher(Y,Fs,window_size,hipass,lopass,numboot)
% Takes the time x channel data matrix Y, filters out a single frequency
% band, and computes the power envelope correlation between each channel
% This is done over time by considering the data in discrete windows
% Input is: 
%           Y  - data matrix, time points by channels
%           Fs - Sampling Frequency in Hz
%           window_size - size of analysis window in seconds
%           hipass - lower bound of filter
%           lopass - upper bound of filter
%           numboot- number of bootstrapped random samples to compare (make zero if no bootstrap wanted for testing)
% Output is:
%            all_coher (the coherence matrix) and           
%            If numboot > 0 then these are expressed in monteCarlo p-value

%--------------------------------------------------------------------------
%% Initialize some variables
n_chan = size(Y,2);  %number of channels
step = Fs*window_size; %Size of each coherence window to analyze in seconds;
n_steps = floor(size(Y,1)/step); %Count how many windows overall
all_coher = zeros(n_steps,n_chan,n_chan);  %Allocate memory for output

fprintf(['The data is ' num2str(size(Y,1)/Fs) ' seconds long. ' '\n']);
fprintf(['Window size is ' num2str(window_size) ' seconds long. ' '\n']);



%--------------------------------------------------------------------------
%% Filter to get the frequency of intereat
order = 3;
bandpass = [(hipass*2)/Fs, (lopass*2)/Fs];          % Bandwidth of bandpass filter
[Bbp,Abp] = butter(order,bandpass);                 % Generation of Xth order Butterworth highpass filter
for c = 1:size(Y,2)
    Y(:,c) = filtfilt(Bbp,Abp,(Y(:,c)));       % Butterworth bandpass filtering of YY
end



%-------------------------------------------------------------------------
%% Compute the correlation between channels 
Y_high_pow = abs(hilbert(Y)).^2;                   % Take the instantaneous angle of the filtered time series in radians using a hilbert transform

frame_count = 0;  %Counter for each time step
for i_time = 1:step:size(Y,1)-step+1                  % for each time window
    tic                                             % start a clock
    range = i_time:i_time+step-1;                     % find the time range based on the current window
    frame_count = frame_count + 1;                  % increment counter
    fprintf(['Processing window ' num2str(frame_count) ' of ' num2str(n_steps) ' windows.  ']);
    
    if numboot == 0    
        all_coher(frame_count,:,:) = corr(Y_high_pow(range,:));  
    else
        
        for i_P1 = 1:n_chan                             % cycle through each channel and take out the data
            seg_P1 = Y_high_pow(range,i_P1);          % data segment
            for i_P2 = 1:n_chan                         % cycle for the channel with which to compare  
                if i_P2 == i_P1                         % in order to avoid rounding errors 
                    all_coher(frame_count,i_P1,i_P2) = 1;
                elseif i_P2 < i_P1                      % Compute the phase coherence and record the values
                    seg_P2 = Y_high_pow(range,i_P2);  % data segment
                    obs_corr = corr(seg_P1,seg_P2);
                    boot_starts = randi([1 size(Y,1)-step],numboot,1);
                    corr_boots = zeros(numboot,1);
                    for i_boot = 1:numboot                                          %for each bootstrap pick a random time point to start the second signal
                        boot_range = boot_starts(i_boot):boot_starts(i_boot)+step-1;  %pick a random start from the list
                        seg_P2_boot = Y_high_pow(boot_range,i_P2);                %grab the data
                        corr_boots(i_boot) = corr(seg_P1,seg_P2_boot);                   %compute the coherence for this boot sample and store in list
                    end
                    all_coher(frame_count,i_P1,i_P2) = length(find(abs(corr_boots) > obs_corr))/numboot;
                end
            end
        end
        
    end
    
    toc                                             %stop clock and display elapsed time
end

