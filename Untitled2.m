%PART-B
%NEURAL NETWORK TRAINING AND TESTING

clear all
close all
clc
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
for qq=1105;   
yy=int2str(qq);
 load(yy,'-mat');
net=newff(inputs_tr,targets,[15 20 15]); %creating a feed forward backprop network with 3 hidden layer of 12,20 &8 resp neuron 
net.trainFcn = 'trainlm';
 
net=train(net,inputs_tr,targets);
 
outputs=net(inputs_tst);
smoothsig = sgolayfilt(outputs,7,21);
outputs=smoothsig;
figure
plot(ecg_data_base_tst,'r');hold on;
plot(outputs,'g')
title('red_is_ECG_test,green_is_PREDICTED_ECG')
for n=0:(length(outputs)/200)-1
figure,
plot(1:200,ecg_data_base_tst(200*n+1:200*n+200),1:200,2*outputs(200*n+1:200*n+200));
end
end
