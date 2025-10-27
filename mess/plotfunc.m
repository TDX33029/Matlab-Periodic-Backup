close all;
baseInfo = csvread('basementInfo.csv');
NumberOfBasement = baseInfo(1,1);
data_original = csvread('OriginalLocation.csv');
data_processed1 = csvread('TDOA_PROCESSED.csv');
data_processed2 = csvread('AOA_PROCESSED.csv');
lengthx = size(data_original);

figure(1);
hold on;
grid on;
axis equal;
set(gca,'ylim',[0 50], 'FontSize', 16);
set(gca,'xlim',[0 100], 'FontSize', 16);
datag = [];
datab = [];
cir = [];
% 
for cir0 = 1:lengthx
    cir = [cir;cir0];
    datag = [data;sqrt((data_original(cir0,1)-data_processed1(cir0,1))^2+(data_original(cir0,2)-data_processed1(cir0,2))^2)];
    datab = [data;sqrt((data_original(cir0,1)-data_processed1(cir0,1))^2+(data_original(cir0,2)-data_processed1(cir0,2))^2)];
    %data = [data;(sqrt((data_original(cir0,1)-data_processed1(cir0,1))^2+(data_original(cir0,2)-data_processed1(cir0,2))^2)+sqrt((data_original(cir0,1)-data_processed2(cir0,1))^2+(data_original(cir0,2)-data_processed2(cir0,2))^2))/2];
end

plot(cir(:,1),data(:,1),'--','LineWidth', 1);
xlabel('Sample Point', 'FontSize', 16);
ylabel('Location Error (cm)', 'FontSize', 16);

