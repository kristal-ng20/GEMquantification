%% Plot to compare cdf plots

fig = figure();
ax = axes(fig);
hold on


fh1 = open('HEK_GEM_Hypertonic_5min006_TracksCDFplot.fig');
axesHandlesToAllLines1 = findobj(fh1, 'Type', 'line');
set(axesHandlesToAllLines1(1),'Color','r')
data1=get(gca,'Children'); %get the handle of the line object
data1XData=get(data1,'XData'); %get the x data
data1YData=get(data1,'YData'); %get the y data
data1Data=[data1XData' data1YData']; %join the x and y data on one array nx2

fh2 = open('HEK_GEM_Hypotonic_10min006_TracksCDFplot.fig');
axesHandlesToAllLines2 = findobj(fh2, 'Type', 'line');
set(axesHandlesToAllLines2(1),'Color','b')
data2=get(gca,'Children'); %get the handle of the line object
data2Data=get(data2,'XData'); %get the x data
data2Data=get(data2,'YData'); %get the y data
data2Data=[data2Data' data2Data']; %join the x and y data on one array nx2

fh3 = open('HEK_GEM_sorbitol_30min006_TracksCDFplot.fig');
axesHandlesToAllLines3 = findobj(fh3, 'Type', 'line');
set(axesHandlesToAllLines3(1),'Color','g')
data3=get(gca,'Children'); %get the handle of the line object
data3XData=get(data3,'XData'); %get the x data
data3YData=get(data3,'YData'); %get the y data
data3Data=[data3XData' data3YData']; %join the x and y data on one array nx2

fh4 = open('HEK_GEM_control003_TracksCDFplot.fig');
axesHandlesToAllLines4 = findobj(fh4, 'Type', 'line');
set(axesHandlesToAllLines4(1),'Color','y')
data4=get(gca,'Children'); %get the handle of the line object
data4XData=get(data4,'XData'); %get the x data
data4YData=get(data4,'YData'); %get the y data
data4Data=[data4XData' data4YData']; %join the x and y data on one array nx2

copyobj(fh1.Children.Children,ax);
copyobj(fh2.Children.Children,ax);
copyobj(fh3.Children.Children,ax);
copyobj(fh4.Children.Children,ax);

close(fh1);
close(fh2);
close(fh3);
close(fh4);


xlabel('Deff (um2.s-1)')
ylabel('CDF(Deff)')
title('CDF Plot')
legend('Hypertonic','Hypotonic','Sorbitol','Control','Location','best')


[h,p] = kstest2([data4Data(:)],[data3Data(:)],'Alpha',0.01)
[h,p] = kstest2([data4Data(:)],[data2Data(:)],'Alpha',0.01)
[h,p] = kstest2([data4Data(:)],[data1Data(:)],'Alpha',0.01)
[h,p] = kstest2([data1Data(:)],[data2Data(:)],'Alpha',0.01)
[h,p] = kstest2([data1Data(:)],[data3Data(:)],'Alpha',0.01)
[h,p] = kstest2([data3Data(:)],[data2Data(:)],'Alpha',0.01)



