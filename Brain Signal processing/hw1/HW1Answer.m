function HW1Answer()
%% This is a matlab syntac that reports the answers to Homework 1 for the ELE573 Brain Signal Processing and App for Fall 2023 C Behtom Adeli
% help goes here

load("sampleEEGdata.mat");

%% a) Topoplot of Averaged 0 to 800ms

start_point = find(EEG.times>-20,1); 
end_point = find(EEG.times>840,1);

split_time=(EEG.times(1,end)-EEG.times(1,1))/length(EEG.times);
step_20ms= (20/round(split_time));
step_100ms=(100/round(split_time));

averaged_splits=mean(EEG.data(:,start_point:end_point,:),3);

figure();
id= 1;

for step=step_20ms+1:step_100ms:length(averaged_splits)
    meaned_data =double(mean(averaged_splits(:,step-step_20ms:step+step_20ms),2));
    subplot(3,3,id);
    topoplot(meaned_data,'eloc64C2.txt','EEG','ColorMap','Jet')
    xlabel('head');                                                                                               
    title(['Topoplot of ',num2str((id-1)*100), 'ms']);         
    id=id+1;
end

%% b)

figure();
channels = cell(1,EEG.nbchan);
max_peak = zeros(2,EEG.nbchan);

for chnl=1:size(averaged_splits,1)
    channels{chnl}=averaged_splits(chnl,:);
    [max_peak(1,chnl), max_peak(2,chnl)] = max(channels{chnl});
    max_peak(2,chnl)=EEG.times(max_peak(2,chnl)+start_point);
    subplot(8,8,chnl);
    plot(EEG.times(1,start_point:end_point),channels{chnl})                                                                                            
    title(['Channel ',num2str(chnl)]);

end



figure();
topoplot(double(max_peak(2,:)),'eloc64C2.txt','EEG','ColorMap','Jet');
title(['Topoplot of ARP Max Peak Latencies']);
hcb=colorbar('southoutside');
hcb.Title.String = "ARP Peak Time Latencies in miliseconds";


%% C)

eloc64_fileTable= readtable('eloc64C2.txt');
theta =  table2array(eloc64_fileTable(:,2));
radius =  table2array(eloc64_fileTable(:,3));
[ X , Y ] = pol2cart( deg2rad(theta) , radius );

%{
% get XYZ coordinates in convenient variables
X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];
%}

% initialize distance matrices
eucdist = zeros(EEG.nbchan,EEG.nbchan);



for chnl2=1:EEG.nbchan
    for chnl=1:EEG.nbchan
        eucdist(chnl2,chnl) = sqrt( (X(chnl)-X(chnl2))^2 + (Y(chnl)-Y(chnl2))^2 );
        %{
        eucdist(chnl2,chnl) = sqrt( (EEG.chanlocs(chnl).X-EEG.chanlocs(chnl2).X)^2 + ...
                                            (EEG.chanlocs(chnl).Y-EEG.chanlocs(chnl2).Y)^2 ); ...
                                              (EEG.chanlocs(chnl).Z-EEG.chanlocs(chnl2).Z)^2 );
        %}
    end
end

lo_width = 0.18;
hi_width = 0.28;
eucdist_filtered = zeros(64,64);
figure();
id =1;
for chn1=1:EEG.nbchan
    for chn2=1:EEG.nbchan
        if (eucdist(chn1,chn2)>lo_width && eucdist(chn1,chn2)<hi_width)
            eucdist_filtered(chn1,chn2) = eucdist(chn1,chn2);
        end
    end
    subplot(8,8,id);
    topoplot(double(eucdist_filtered(chn1,:)),'eloc64C2.txt','EEG','ColorMap','Jet');
    id = id +1 ;
    title(['channel' eloc64_fileTable{chn1,4}{1,1}]);
end

    
for chn1=1:EEG.nbchan
    for chn2=1:EEG.nbchan
        if eucdist_filtered(chn1,chn2)
        weigths= (1/eucdist_filtered(chn1,chn2))/sum(eucdist_filtered(chn1,:));
        end
    end
end



disp('it is done.')