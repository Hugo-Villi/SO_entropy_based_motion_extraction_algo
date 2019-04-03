function fig=plot_scatter3D(all_led_info,data_to_plot,number_of_bars)
%creation of the figure object
fig=figure('position',[50,50,1400,700]...
    ,'name','optogait data visualization'); %Set the position and dimension as well as the name

%creation of the plot
%p=plot(all_led_info(1,:))

%creation of a text display
txh=uicontrol(fig,'Style','text','position',[50,650,100,20]); %will display the frame number

%creation of axes
axes('units', 'pixels', ...
    'position', [100 100 1200 530]);    %set position and dimension of the axis (plot)

%Creation of the object Uicontrol slider
uicontrol(fig,'Style','slider'...   %define the type of controller
    ,'position',[25,25,1350,30],... %define the position and dimensions
    'min',0,'max',size(all_led_info,1),...  %define the max/min limits
    'SliderStep',[1/size(all_led_info,1) (1/size(all_led_info,1))*10],...   %define the steps for a click on the arrows and for a click directly on the slider
    'callback',@update);  %callback to the update function

    function update(h,~)    %declare the update function
        test_value=get(h,'Value');  %get the value of the slider
        test_value=round(test_value);   %round the value to have an integer in order to use it to get data from the array
        text_to_disp=data_to_plot(test_value,1);    %get the frame number
        set(h,'Value',test_value);  %reset the slider to the rounded value
        set(txh,'string',num2str(text_to_disp));    %set the text box to display the frame number
        set(txh,'FontSize',12);
        plot(all_led_info(test_value,:)); %actually plot the state of the leds associate to the frame
        ax=gca;  %get the index of the axis to modify them
        ax.XLim=[0 number_of_bars*96];   %set the limit in function of the number of bars of the setup
        %ax.XTick = 0:96:(number_of_bars*96);    %set the tick to respect the disposition of the setup
        ax.XTick = 0:10:(number_of_bars*96);    %more precise
        ax.XMinorTick = 'on' ;  %set minor ticks between the majors one defined at the previous line
        ax.XLabel.String='Leds';   %axis X label
        ax.YLim=[0 1.5];    %Y limits
        ax.YTick=[0,1]; %Y ticks
        ax.YLabel.String='0=not cut, 1=cut';    %axis Y label
    end
end