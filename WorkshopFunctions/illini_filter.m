function [Y_filt] = illini_filter(Y,Fs,hipass,lopass,order,type)
%% illini_filter.m - Created University of Illinois, based on kyle_filter by Kyle Mathewson with edits by Jamie Norton

% [Y_filt] = illini_filter(Y,Fs,hipass,lopass,order,type)
% Takes the time x channel data matrix Y, filters out a single frequency band,
% Input is: 
%           Y  - data matrix, time points by channels
%           Fs - Sampling Frequency in Hz
%           hipass - lower bound of filter
%           lopass - upper bound of filter
%           order - filter order, 
%                      higher orders have steeper roll offs, but
%                      can also be unstable
%           type - 'low' (low pass), 'high' (high pass), or 'band' (bandpass)
%
% Output is:
%            Y_filt - same as Y but not filtered

%--------------------------------------------------------------------------
%% Initialize some variables
Y_filt = Y;  %Allocate memory for output
bandpass = [];
if strcmp(type,'band') == 1
    bandpass = [(hipass*2)/Fs, (lopass*2)/Fs]; 
    fprintf(['Filtering data between ' num2str(hipass) ' Hz and ' num2str(lopass) ' Hz. ' '\n']);
end  
if strcmp(type,'high') == 1
    bandpass = (hipass*2)/Fs;
    fprintf(['Filter removing below ' num2str(hipass) ' Hz' '\n']);
end
if strcmp(type,'low') == 1
    bandpass = (lopass*2)/Fs;
    fprintf(['Filter removing above ' num2str(lopass) ' Hz' '\n']);
end    

%--------------------------------------------------------------------------
%% Filter with a bandpass butterworth
[Bbp,Abp] = butter(order,bandpass);                 % Generation of Xth order Butterworth highpass filter
for c = 1:size(Y,2)
    Y_filt(:,c) = filtfilt(Bbp,Abp,(Y(:,c)));       % Butterworth bandpass filtering of YY
end

%% plot the data
t_data = 1/Fs:1/Fs:(1/Fs)*length(Y);
figure; subplot(2,1,1); plot(t_data,Y); title('Original Data'); ylabel('Voltage (uV)'); xlabel('Time (ms)'); axis tight;
subplot(2,1,2); plot(t_data,Y_filt); title('Filtered Data'); ylabel('Voltage (uV)'); xlabel('Time (ms)'); axis tight;