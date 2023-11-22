%PART –A 
%CREATING DATABASE

close all
clear all
clc
rr=1099;
for q=100:130;  %no. of data from different patient
y=int2str(q);
load(y);
rr=rr+1;  %for saving database serially for all the available person
zz=int2str(rr);
close all;
nos_of_signals=length(val(:,1));
v=1;                      % "1" for ecg; "2" for bp; "3" for eeg...taking as base signal
for n=1:nos_of_signals;
sig_data=val(n,:);
t = 1:length(sig_data);
%sig_data= (sig_data/max(abs (sig_data))); %normalization of signal
smoothsig = sgolayfilt(sig_data,7,21);  % smoothing the signal using sgolayfit
figure
plot(t(1:1000),sig_data(1:1000),'r',t(1:1000),smoothsig(1:1000),'b'); grid on
axis tight;
xlabel('Samples'); ylabel('Voltage(mV)');
legend('Noisy  Signal','sgolay Filtered Signal')
title('Filtering Noisy  Signal')

vall(n,:)=smoothsig;
t(n,:)=1:length(vall(n,:));
%normalizing the  filtered data signal
max_amp=max(abs(vall(n,:)));
normlized_sig(n,:)=(vall(n,:))/max_amp;
figure
plot(t(n,1:1000),normlized_sig(n,1:1000),'b'); %ploting the  normalized raw signal
title('filtered normalized signal')
grid on

% to performing 5pt differentiation
k(n,:)=[0,0,normlized_sig(n,:),0,0];
for i=1:1:numel(k(n,:))-4
diff5pt(n,i)=((-k(n,(i+4))+(8*((k(n,(i+3)))))-(8*k(n,(i+1)))+k(n,(i)))*1000)/12;
end
figure
stem(t(n,1:1000),diff5pt(n,1:1000),'r'); hold on;%ploting the 5point diff  signal
title('5point diff  signal')
grid on
%taking absolute of the 5pt diff signal to get more exact peak of signal
abs_d(n,:)=abs(diff5pt(n,:));
stem(t(n,1:1000),abs_d(n,1:1000),'b');%ploting the  absolute 5point diff  signal
title(' absolute 5point diff  Signal');
grid on;


%to find the highest peak value in the signal,not considering the peak occurs at 1st,2nd ,2nd_last & last sample
hig_peak(n,:)=max(abs_d(n,3:length(abs_d(n,:))-2));  % higest peak value
%to find the time of this highest peak value in the signal
for i=1:1:length(abs_d(n,:))
if(hig_peak(n,:) == abs_d(n,i))
hig_peak_time(n,:)= i;
end
end

%to find time around which there may be a max peak
y1=0; %clearing the array y1
x1=hig_peak_time(n,:)  : 200 : length(vall(n,:)); %x1 goes in increasing order
y= hig_peak_time(n,:)  :-200 : 1 ;               %y goes in decreasing order
for i=0:1:abs(length(y)-2)
if(length(y)==1)
i=0;
else
y1(i+1)=y((length(y)-i));  %making the y as y1 that goes in increading order till second last value of y
end
end
if(length(y)==1)
z=x1;
else
z=[y1,x1]; %the array z is having the value of  approx time in increasing order around which there may be a max-peak
end
for i=1:1:length(z)
time(n,i)=z(i);% matrix 'time' is having the value of approx time in increasing order around which there may be a max-peak
end


%the value of time in increasing order at which there lies a peak
for mask=1:4
for j=1:1:length(time(n,:))
t2=time(n,j);
lsv=(t2-100);%least significant value
msv=(t2+100);%most significant value
if(lsv<=0)
lsv=1;
end
if(msv>length(vall(n,:)))
msv=length(vall(n,:));
end
for i= lsv:1:msv
if(abs_d(n,time(n,j))< abs_d(n,i))
time(n,j)=i;
end
end
end
end

%removing the repeatition peak time
c=0;
for i=1:length(time(n,:))-1
if (abs( time(n,i+1) - time(n,i)  > 100))
c=c+1;
time1(n,c)=time(n,i);
end
end
end

for n=1:nos_of_signals;

%the value of time (of selected 'v'signal) in increasing order at which there lies a peak
for i=1:length(time1(v,:))
if(time1(v,i)>0)
time2(v,i)=time1(v,i);
end
end
time2=time2;    %since the 5point diff lag the peak compare to raw signal by 4 sample.
figure
stem(time2(v,:),'r-','MarkerFaceColor','b')
title(' time at which there is a peaks in the  ecg signal ');


%making the data base of signal taking the time array time at which there
%is a peaks in the  ecg signal in increasing order.
a=2; d_b=0;  m=1; w=1;

for a=2:4:length(time2(v,:))-1
d_b= normlized_sig(n,(time2(v,a-1):time2(v,a+1)));
for l=1:length(d_b)
data_base_1(1,w)=d_b(l);
w=w+1;
end
a=a+1;
if(a==length(time2(v,:)))
break;
end
d_b= normlized_sig(n,(time2(v,a-1):time2(v,a+1)));
for l=1:length(d_b)
data_base(1,m)=d_b(l);
m=m+1;
end
a=a+1;
if(a==length(time2(v,:)))
break
end
d_b= normlized_sig(n,(time2(v,a-1):time2(v,a+1)));
for l=1:length(d_b)
data_base(1,m)=d_b(l);
m=m+1;
end
a=a+1;
if(a==length(time2(v,:)))
break
end
d_b= normlized_sig(n,(time2(v,a-1):time2(v,a+1)));
for l=1:length(d_b)
data_base(1,m)=d_b(l);
m=m+1;
end

end
if( n==1)

ecg_data_base_tr=(data_base);
ecg_data_base_tst=(data_base_1);
elseif(n==2)
bp_data_base_tr=(data_base);
bp_data_base_tst=(data_base_1);
elseif(n==3)
eeg_data_base_tr=(data_base);
eeg_data_base_tst=(data_base_1);
elseif( n==4)
sg4_data_base_tr=(data_base);
sg4_data_base_tst=(data_base_1);
elseif( n==5)
sg5_data_base_tr=(data_base);
sg5_data_base_tst=(data_base_1);
elseif( n==6)
sg6_data_base_tr=(data_base);
sg6_data_base_tst=(data_base_1);
else
sg7_data_base_tr=(data_base);
sg7_data_base_tst=(data_base_1);
end
end



for a=2:1:length(time2(v,:))-(length(time2(v,:))-4)

figure
subplot(2,2,4);
plot(t(1,(time2(v,a-1):time2(v,a+1))),abs_d(1,(time2(v,a-1):time2(v,a+1))),'-r');hold on;
plot(t(n,(time2(v,a-1):time2(v,a+1))),abs_d(n,(time2(v,a-1):time2(v,a+1))),'-b');
title('red is ecg abs_d sig blue is 2nd sig')

subplot(2,2,3);
plot(t(1,(time2(v,a-1):time2(v,a+1))),diff5pt(1,(time2(v,a-1):time2(v,a+1))),'-r');hold on;
plot(t(n,(time2(v,a-1):time2(v,a+1))),diff5pt(n,(time2(v,a-1):time2(v,a+1))),'-b');
%legend('red is diff5pt sig',' blue is 2nd sig')
title('red is ecg diff5pt sig blue is 2nd sig')

subplot(2,2,2);
plot(t(1,(time2(v,a-1):time2(v,a+1))),vall(1,(time2(v,a-1):time2(v,a+1))),'-r');hold on;
plot(t(n,(time2(v,a-1):time2(v,a+1))),vall(n,(time2(v,a-1):time2(v,a+1))),'-b');
%legend('red is filtred ecg sig',' blue is 2nd sig')
title('red is filtred ecg sig blue is 2nd sig')

subplot(2,2,1);
plot(t(1,(time2(v,a-1):time2(v,a+1))),normlized_sig(1,(time2(v,a-1):time2(v,a+1))),'-r');hold on;
plot(t(n,(time2(v,a-1):time2(v,a+1))),normlized_sig(n,(time2(v,a-1):time2(v,a+1))),'-b');
%legend('red is normalized filtred ecg  sig',' blue is 2nd sig')
title('red is ecg normalized filtred ecg  sig blue is 2nd sig')
end

%save the data base of available signals to a file named zz
if( n==1)

inputs_tr=ecg_data_base_tr; %training data
inputs_tst=ecg_data_base_tst;%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','ecg_data_base_tst','time2','inputs_tr','inputs_tst','targets') ;

elseif(n==2)
inputs_tr=bp_data_base_tr; %training data
inputs_tst=bp_data_base_tst;%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','bp_data_base_tr','ecg_data_base_tst','bp_data_base_tst','time2','inputs_tr','inputs_tst','targets') ;

elseif(n==3)

inputs_tr=[bp_data_base_tr;eeg_data_base_tr]; %training data
inputs_tst=[bp_data_base_tst;eeg_data_base_tst];%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','bp_data_base_tr','eeg_data_base_tr','ecg_data_base_tst','bp_data_base_tst','eeg_data_base_tst','time2','inputs_tr','inputs_tst','targets') ;

elseif( n==4)
inputs_tr=[bp_data_base_tr;eeg_data_base_tr;sg4_data_base_tr]; %training data
inputs_tst=[bp_data_base_tst;eeg_data_base_tst;sg4_data_base_tst];%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','bp_data_base_tr','eeg_data_base_tr','sg4_data_base_tr','ecg_data_base_tst','bp_data_base_tst','eeg_data_base_tst','sg4_data_base_tst','time2','inputs_tr','inputs_tst','targets') ; % function form to save variable of workspace

elseif( n==5)
inputs_tr=[bp_data_base_tr;eeg_data_base_tr;sg4_data_base_tr;sg5_data_base_tr]; %training data
inputs_tst=[bp_data_base_tst;eeg_data_base_tst;sg4_data_base_tst;sg5_data_base_tst];%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','bp_data_base_tr','eeg_data_base_tr','sg4_data_base_tr','sg5_data_base_tr','ecg_data_base_tst','bp_data_base_tst','eeg_data_base_tst','sg4_data_base_tst','sg5_data_base_tst','time2','inputs_tr','inputs_tst','targets') ; % function form to save variable of workspac

elseif( n==6)
inputs_tr=[bp_data_base_tr;eeg_data_base_tr;sg4_data_base_tr;sg5_data_base_tr;sg6_data_base_tr]; %training data
inputs_tst=[bp_data_base_tst;eeg_data_base_tst;sg4_data_base_tst;sg5_data_base_tst;sg6_data_base_tst];%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','bp_data_base_tr','eeg_data_base_tr','sg4_data_base_tr','sg5_data_base_tr','sg6_data_base_tr','ecg_data_base_tst','bp_data_base_tst','eeg_data_base_tst','sg4_data_base_tst','sg5_data_base_tst','sg6_data_base_tst','time2','inputs_tr','inputs_tst','targets') ; % function form to save variable of workspace

else
inputs_tr=[bp_data_base_tr;eeg_data_base_tr;sg4_data_base_tr;sg5_data_base_tr;sg6_data_base_tr;sg7_data_base_tr]; %training data
inputs_tst=[bp_data_base_tst;eeg_data_base_tst;sg4_data_base_tst;sg5_data_base_tst;sg6_data_base_tst;sg7_data_base_tst];%test data
targets=ecg_data_base_tr;
save(zz,'-mat','ecg_data_base_tr','bp_data_base_tr','eeg_data_base_tr','sg4_data_base_tr','sg5_data_base_tr','sg6_data_base_tr','sg7_data_base_tr','ecg_data_base_tst','bp_data_base_tst','eeg_data_base_tst','sg4_data_base_tst','sg5_data_base_tst','sg6_data_base_tst','sg7_data_base_tst','time2','inputs_tr','inputs_tst','targets') ; % function form to save variable of workspace
end
end


%PART-B
%NEURAL NETWORK TRAINING AND TESTING

%clear all
%close all
%clc
% for qq=1105;   
% yy=int2str(qq);
%  load(yy,'-mat');
% net=newff(inputs_tr,targets,[15 20 15]); %creating a feed forward backprop network with 3 hidden layer of 12,20 &8 resp neuron 
% net.trainFcn = 'trainlm';
%  
% net=train(net,inputs_tr,targets);
%  
% outputs=net(inputs_tst);
% smoothsig = sgolayfilt(outputs,7,21);
% outputs=smoothsig;
% figure
% plot(ecg_data_base_tst,'r');hold on;
% plot(outputs,'g')
% title('red_is_ECG_test,green_is_PREDICTED_ECG')
% for n=0:(length(outputs)/200)-1
% figure,
% plot(1:200,ecg_data_base_tst(200*n+1:200*n+200),1:200,2*outputs(200*n+1:200*n+200));
% end
% end
