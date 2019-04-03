clearvars -except files k marker_set SAVE_DATA save_file_name first_frame g
file_name = uigetfile('*.c3d');
acq = btkReadAcquisition(file_name); %read the file
markers = btkGetMarkers(acq);   %get the information of the markers in a structure
labels = fieldnames(markers);   %get the labels of the markers
markers_values = btkGetMarkersValues(acq);
for i=1:5
    order=4;
    f=btkGetPointFrequency(acq);
    cutoff_freq=i*2;
    norm_cutoff_freq=cutoff_freq/(f/2);
    [b,a]=butter(order,norm_cutoff_freq,'low');
    z=markers_values(:,3);
    subplot(3,2,1);
    plot(z)
    subplot(3,2,i+1);
    plot(filtfilt(b,a,z));
    title(i*2);
end
