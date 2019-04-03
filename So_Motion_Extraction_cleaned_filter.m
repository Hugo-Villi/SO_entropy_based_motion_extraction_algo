%%acquisition of the needed data from the motion capture.
clearvars
myFolder = 'D:\stage\Sofamehack2019\test';  %define the folder where the files to compute are
files=get_file_names_c3d(myFolder); %Get the names of the files inside the folder
for z=1:2
    for k=1:size(files,1)
        %clearvars -except files k marker_set SAVE_DATA save_file_name first_frame g z %clear all the variables except the one needed toexecute the loop over all the files
        file_name = files(k).name;  %get the name of the file that will be treated in this iteration
        acq = btkReadAcquisition(file_name); %read the file
        labels = fieldnames(btkGetMarkers(acq));   %get the labels of the markers
        markers_values = btkGetMarkersValues(acq); %give an array filled with the coordinates of the marker (x,y,z of first marker, x,y,z of second, and so on)
        events=btkGetEvents(acq);   %get the events stored in the c3d file (rigth foot strike...)
        first_frame(k)=btkGetFirstFrame(acq);   %Get the number of the first frame to correct the delay between the events detected by the algo and the event directly in the files
        %{ 
        %filter for the test files generated with a basketball and a sandbag
        order=4;    %order of the filter
        f=btkGetPointFrequency(acq);    %get the frequency of acquisition of the file
        cutoff_freq=10; %set the desired cut off frequency
        norm_cutoff_freq=cutoff_freq/(f/2); %give the normalized cutoff frequency to input in the butter function
        [b,a]=butter(order,norm_cutoff_freq,'low'); %gives the transfer function coefficients for the lowpass digital butterworth filet
        if z==2;    %makes a test only to compare filtered and unfiltered results
            for i=1:size(markers_values,2)
                markers_values(:,i)=filtfilt(b,a,markers_values(:,i)); %applies  a zero-phase delay filtering to each x,y,z component of each marker
            end
        end
        subplot(2,1,1)  %create a subplot, where multiple figure can be drawn in one window
        plot(markers_values(:,3))   %plot the z components of the first marker
        hold on     %the hold on will allow to draw mutliple line on the same figure
        %}
        %%calculation of the displacement vector
        num_row = size(markers_values,1);   %a row = a frame
        num_col = size(markers_values,2);   %this is equal to 3*number of markers
        displacement=zeros(num_row,num_col);    %preallocation
        for i=1:(num_row-1) %-1 otherwise it will crash trying to reach the last frame + 1
            for j=1:num_col
                displacement(i,j)=markers_values(i+1,j)-markers_values(i,j);    %simple soustraction between frame+1-frame
            end
        end
        %give only displacements for x, y and z, to ease the next steps (max/min
        %detection and discretization)
        for i=3:3:num_col
            displacement_x(:,i/3)=displacement(:,i-2);
        end
        for i=3:3:num_col
            displacement_y(:,i/3)=displacement(:,i-1);
        end
        for i=3:3:num_col
            displacement_z(:,i/3)=displacement(:,i);
        end
        %%displacement histogram
        n = 500; %number of discretizations levels, this setting may have importance on the results
        maximum=max(displacement);  %gives the max values for all the column hence for all the x,y,z coordinates respectively
        minimum=min(displacement);  %same for min
        %the following part get the maximum and minimums of x,y and z displacement
        %to set the range of the histograms
        max_x = max(displacement_x,[],'all');
        min_x = min(displacement_x,[],'all');
        
        max_y = max(displacement_y,[],'all');
        min_y = min(displacement_y,[],'all');
        
        max_z = max(displacement_z,[],'all');
        min_z = min(displacement_z,[],'all');
        
        %creating the histograms
        discretization_x=linspace(min_x,max_x,n);   %the function linspace creates a vector of values envenly distributed along a range
        for i=1:size(displacement_x,1)    %will generate an histogram for each frame for the x coordinates
            histogram_x(i,:)=histcounts(displacement_x(i,:),discretization_x);
        end
        discretization_y=linspace(min_y,max_y,n);
        for i=1:size(displacement_y,1)
            histogram_y(i,:)=histcounts(displacement_y(i,:),discretization_y);
        end
        discretization_z=linspace(min_z,max_z,n);
        for i=1:size(displacement_z,1)
            histogram_z(i,:)=histcounts(displacement_z(i,:),discretization_z);
        end
        
        
        %computing the probabilities
        for i=1:size(histogram_x)
            for j=1:size(histogram_x,2)
                prob_x(i,j)=histogram_x(i,j)/(num_col/3);
            end
        end
        for i=1:size(histogram_x)
            for j=1:size(histogram_x,2)
                prob_y(i,j)=histogram_y(i,j)/(num_col/3);
            end
        end
        for i=1:size(histogram_x)
            for j=1:size(histogram_x,2)
                prob_z(i,j)=histogram_z(i,j)/(num_col/3);
            end
        end
        
        
        %filling the C(r,s) matrix
        %with the same global range as the histograms
        %create the matrices for each marker
        [matrix_r_s_global_x,matrix_r_s_x,Cx,proba_r_s_x]=compute_C(histogram_x,displacement_x,discretization_x,labels);
        [matrix_r_s_global_y,matrix_r_s_y,Cy,proba_r_s_y]=compute_C(histogram_y,displacement_y,discretization_y,labels);
        [matrix_r_s_global_z,matrix_r_s_z,Cz,proba_r_s_z]=compute_C(histogram_z,displacement_z,discretization_z,labels);
        
        %compute the mutual information for each dimension and sum it to get the
        %total mutual information I
        Ix=mutual_info(displacement,Cx,prob_x);
        Iy=mutual_info(displacement,Cy,prob_y);
        Iz=mutual_info(displacement,Cz,prob_z);
        I=Ix+Iy+Iz;
        
        file_name=erase(file_name,'.c3d');
        save_file_name(k)=string(file_name);
        SAVE_DATA.(file_name).Ix=Ix;
        SAVE_DATA.(file_name).Iy=Iy;
        SAVE_DATA.(file_name).Iz=Iz;
        SAVE_DATA.(file_name).I=I;
        %if isa(events,'struct')==1
            %events=events.General_impact;
        %end
        %SAVE_DATA.(file_name).forceplate=events;
        subplot(2,1,2)
        plot(Iz)
        hold on
    end
end
hold off
%legend('raw','filtered')
%{
%set=["all","legs","shank","foot"];
window_size=20;
threshold=0;
clear comparison
comparison(1,1)=string('window size');
comparison(1,2)=window_size;
comparison(1,3)=string('threshold');
comparison(1,4)=threshold;
b=2;
for j=1:size(save_file_name,2)
    clear keyposes_temp keyposes_frames I_localized I_trimmed
    name=char(save_file_name(1,j));
    I=SAVE_DATA.(name).I;
    threshold=0;
    [keyposes_temp,I_localized,I_trimmed]=keyposes_detection(I,window_size,threshold);
    while size(keyposes_temp,2)>3
        [keyposes_temp,I_localized,I_trimmed]=keyposes_detection(I,window_size,threshold);
        threshold=threshold-1;
    end
    
    keyposes_frames(1:size(keyposes_temp,2))=keyposes_temp;
    SAVE_RESULTS.(name).keyposes_frames=keyposes_frames;
    SAVE_RESULTS.(name).I_localized=I_localized;
    SAVE_RESULTS.(name).I_trimmed=I_trimmed;
    %get the timing for the i_trimmed
    frames_global=keyposes_frames+first_frame(j);
    time_cine=frames_global/200;
    SAVE_RESULTS.(name).time_cine=time_cine;
    comparison(b,1)=string(name);
    b=b+1;
    comparison(b,1)='events: forceplate';
    b=b+1;
    comparison(b,1:size(SAVE_DATA.(name).forceplate,2))=SAVE_DATA.(name).forceplate;
    b=b+1;
    comparison(b,1)=('events: algo');
    b=b+1;
    comparison(b,1:size(time_cine,2))=time_cine;
    b=b+2;
    time_cine=[];
end
comparison=fillmissing(comparison,'constant',string(' '));


%plot(I_trimmed);
%scatter_3D_plot(markers_values,keyposes_frames)
%legend()

%}
%[frame_contact,frame_off,frame_contact2,frame_off2]=detect_event_forceplate(acq);
