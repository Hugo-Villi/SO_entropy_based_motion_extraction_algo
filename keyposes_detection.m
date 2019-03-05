function [keyposes_frames,I_localized] = keyposes_detection(I,window_size,threshold)
%Keypose detection
%localization of mutual information to emphasize changes
I_sliding_mean=movmean(I,window_size);  %total window size, will take the mean of (window_size) frames
for i=1:size(I,2)
    I_localized(i)=I(i)/I_sliding_mean(i);
end
%plot(I_localized,'DisplayName',num2str(window_size));
%legend()
%hold on
%pick low minima values
%with a threshold trimming values horizontally
for i=1:size(I_localized,2)
    if I_localized(i)>threshold
        I_trimmed(i)=threshold; %if the value is higher than the threshold, it is replaced by the threshold itself
    else
        I_trimmed(i)=I_localized(i);    %otherwise it is not changed
    end
    
end

%identify local through
j=1;
for i=1:size(I_trimmed,2)
    if I_trimmed(i)==threshold  %if the value is equal to the threshold it means that this value is not in a peak
        negative_peaks(j)=1;    %a value equal to 1 means there the value i is not in a peak
        j=j+1;
    else
        negative_peaks(j)=0;
        j=j+1;
    end
end
negative_peaks(1)=1;    %set a 1 at the beginning and at the end to avoid errors is the further part of the code
negative_peaks(end)=1;

%find the start and ends of the peak
start=1;
stop=1;
k=1;
i=1;
j=1;
while i<size(I_trimmed,2)-1 %to go through the whole array
    while negative_peaks(i)~=0  %will look for a 0. Finding a 0 mean that the start of the peak have benn reached
        i=i+1;
        if i==size(I_trimmed,2) %if the end of the array is reached the while loops are exited
            break
        end
    end
    start(k)=i-1;   %the start frame is stored
    j=start(k)+1;   %j is the starting point of the loop that will look for the end of the peak
    while negative_peaks(j)~=1  %will look for a 1. Finding a 1 mean that the end of the peak have benn reached
        j=j+1;
    end
    stop(k)=j;  %the end frame is stored
    i=stop(k);  %set i as the next starting point for the loop that will search for the e4next start
    k=k+1;  %to increment start/stop arrays
end

%Get the min of the peak
for i=1:size(start,2)
    peaks_test(i)=min(I_trimmed(1,start(1,i):stop(1,i)));   %use the start and stop info to look for min values of the peak
    keyposes_frames(i)=find(I_trimmed==peaks_test(i));  %use of the 'find' function to get the frame where the minimal occur
end
end