function scatter_plot = scatter_3D_plot(markers_values,keyposes_frames)
fig=figure('position',[50,50,1200,700]...
    ,'name','optogait data visualization'); %Set the position and dimension as well as the name

uicontrol(fig,'Style','slider'...   %define the type of controller
    ,'position',[25,25,1150,30],... %define the position and dimensions
    'min',0,'max',size(keyposes_frames,2),...  %define the max/min limits
    'SliderStep',[1/size(keyposes_frames,2) (1/size(keyposes_frames,2))*10],...   %define the steps for a click on the arrows and for a click directly on the slider
    'callback',@update);

    function update(h,~)    %declare the update function
        k=1;
        for i=1:3:size(markers_values,2)
            x(:,k)=markers_values(:,i);
            y(:,k)=markers_values(:,i+1);
            z(:,k)=markers_values(:,i+2);
            k=k+1;
        end
        
        for i=1:12
            for k=1:(size(x,2))
                size_of_points(k,i)=i;
                Color(i,1)=i*0.08;
                Color(i,2)=i*0.08;
                Color(i,3)=i*0.08;
            end
        end
        test_value=get(h,'Value');  %get the value of the slider
        test_value=round(test_value);   %round the value to have an integer in order to use it to get data from the array
        %text_to_disp=data_to_plot(test_value,1);    %get the frame number
        set(h,'Value',test_value);  %reset the slider to the rounded value
        %set(txh,'string',num2str(text_to_disp));    %set the text box to display the frame number
        %set(txh,'FontSize',12);
        l=2;
        k=1;
        color_event=[1,1,0];
        frame=keyposes_frames(test_value)
        if frame-100<1
            start=1
        else
            start=frame-100
        end
        if frame+100>size(x,1)
            stop=size(x,1)
        else
            stop=frame+100
        end
        for around=start:20:stop
            l=l+1;
            k=k+1;
            temp_color=Color(k,:)
            for i=1:size(x,2)
                
                if around==frame
                    scatter_plot=scatter3(x(around,i),y(around,i),z(around,i),l*3,color_event,'filled');
                else
                    scatter_plot=scatter3(x(around,i),y(around,i),z(around,i),l*3,temp_color,'filled');
                    axis('equal');
                    %view(2)
                end
                hold on
            end
        end
        hold off
    end
end