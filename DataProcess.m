close all;
%--  Model 2 --%
disp('Reading Config');
TOAEnable = 1;
AOAEnable = 1;
TDOAEnable = 1;

syms x;
%temp value
data_toa = [];
data_aoa = [];
data_tdoa = [];
data_output = [];
%

list = [];
baseInfo = csvread('basementInfo.csv');
NumberOfBasement = baseInfo(1,1);
basex = csvread('basementInfo.csv',1,0,[1,0,NumberOfBasement,0]);
basey = csvread('basementInfo.csv',1,1,[1,1,NumberOfBasement,1]);
NumberOfBasement = baseInfo(1,1);
if NumberOfBasement == 2 || NumberOfBasement == 1
    disp('Tow Few Basements !');
    return;
end

for NumberOfCirculate1 = 1:NumberOfBasement
    if NumberOfCirculate1 ~= NumberOfBasement-1 && NumberOfCirculate1 ~= NumberOfBasement
        list = [list;NumberOfCirculate1,NumberOfCirculate1+1,NumberOfCirculate1+2];
    elseif NumberOfCirculate1 ~= NumberOfBasement
        list = [list;NumberOfCirculate1,NumberOfCirculate1+1,1];
    else
        list = [list;NumberOfCirculate1,1,2];
    end
end
% disp(list);
figure(2);
grid on;
hold on;
xlabel('Displacement X (m)');
ylabel('Displacement Y (m)');
ylim([0,2000]);
axis equal;
set(gca,'ylim',[0 2000]);
set(gca,'xlim',[-1000 1000]);


if TOAEnable == 1   %TOA SETUP---------------------------------------------------------------------
disp('Printing TOA Data');
    data = csvread('TOA_RES.csv');
    lengthx = size(data);
    for cir0 = 1:NumberOfBasement
        label = sprintf('RU%d',cir0);
        text(basex(cir0)+20, basey(cir0)+20, label, 'FontSize', 8, 'Color', 'red');
    end
    scatter(basex,basey,60,'g');

    for cir1 = 1:NumberOfBasement
        circle(cir1) = rectangle('Position', [basex(cir1)-data(1,cir1),basey(cir1)-data(1,cir1),1,1],...
              'Curvature', [1, 1],...
              'EdgeColor', 'b',...
              'LineWidth', 0.2);
    end


    point = [0,0];

    for cir1 = 1:lengthx
        sum = [0,0];
        for cir2 = 1:NumberOfBasement
            temp = IntersectionCheck(basex(list(cir2,1)),basey(list(cir2,1)),data(cir1,list(cir2,1)),basex(list(cir2,2)),basey(list(cir2,2)),data(cir1,list(cir2,2)));
            if( temp == 0 || temp == 1)
                point(1) = basex(list(cir2,1))+((basex(list(cir2,2))-basex(list(cir2,1)))/(data(cir1,list(cir2,1))+data(cir1,list(cir2,2))))*data(cir1,list(cir2,1));
                point(2) = basey(list(cir2,1))+((basey(list(cir2,2))-basey(list(cir2,1)))/(data(cir1,list(cir2,1))+data(cir1,list(cir2,2))))*data(cir1,list(cir2,1));
            elseif(temp == 2)
                d = DistantGet(basex(list(cir2,1)),basey(list(cir2,1)),basex(list(cir2,2)),basey(list(cir2,2)));
                m = (data(cir1,list(cir2,1))^2+d^2-data(cir1,list(cir2,2))^2)/(2*d);
                h = sqrt(data(cir1,list(cir2,1))^2-m^2);
                point_m(1) = basex(list(cir2,1))+((basex(list(cir2,2))-basex(list(cir2,1)))/d)*m;
                point_m(2) = basey(list(cir2,1))+((basey(list(cir2,2))-basey(list(cir2,1)))/d)*m;
                point_t(1) = point_m(1)+h*(basey(list(cir2,2))-basey(list(cir2,1)))/d;
                point_t(2) = point_m(2)-h*(basex(list(cir2,2))-basex(list(cir2,1)))/d;
                point_t(3) = point_m(1)-h*(basey(list(cir2,2))-basey(list(cir2,1)))/d;
                point_t(4) = point_m(2)+h*(basex(list(cir2,2))-basex(list(cir2,1)))/d;

                d3(1) =  DistantGet(point_t(1),point_t(2),basex(list(cir2,3)),basey(list(cir2,3)));
                d3(2) =  DistantGet(point_t(3),point_t(4),basex(list(cir2,3)),basey(list(cir2,3)));

                if(unsign(data(cir1,list(cir2,3)) - d3(1)) > unsign(data(cir1,list(cir2,3)) - d3(2)))
                point(1) = point_t(3);
                    point(2) = point_t(4);
                else
                    point(1) = point_t(1);
                    point(2) = point_t(2);
                end
            else
                disp('Error');
                return;
            end
%             if imag(point(1)) ~= 0
%                 point(1) = sqrt(real(point(1))^2+imag(point(1))^2);
%             end
%             if imag(point(2)) ~= 0
%                 point(2) = sqrt(real(point(2))^2+imag(point(2))^2);
%             end
            sum(1) = sum(1) +  point(1);
            sum(2) = sum(2) +  point(2);
%             scatter(point(1),point(2),10,'b');     %Print TOA Details
        end
%         disp(sum(1));
%         pause(0.2);
        scatter(sum(1)/NumberOfBasement,sum(2)/NumberOfBasement,10,'r');
        
        data_toa = [data_toa;sum(1)/NumberOfBasement,sum(2)/NumberOfBasement];
  
        %pause(0.01);
    end
    for cir0 = 1:NumberOfBasement
        set(circle(cir0),'Position',[basex(cir0),basey(cir0), 0, 0]);
    end
    csvwrite('TOA_PROCESSED.csv',data_toa);
end %TOA END------------------------------------------------------------------------------------------
disp('Complete');

list = [];
b = [];

for cir0 = 1:NumberOfBasement
   if cir0 ~= NumberOfBasement
       list = [list;cir0,cir0+1];
   else
       list = [list;cir0,1];
   end
end


if AOAEnable == 1 %AOA SETUP--------------------------------------------------------------------------
disp('Printing AOA Data');
    data = csvread('AOA_RES.csv');
    lengthx = size(data);
    for cir1 = 1:lengthx
        for cir0 = 1:NumberOfBasement
            k(cir0) = tan(data(cir1,cir0));
            b(cir0) = basey(cir0) - k(cir0).*basex(cir0);
        end
        for cir0 = 1:NumberOfBasement
            y3 = @(x)k(cir0).*(x-basex(cir0))+basey(cir0);
        end
        sumx = 0;
        sumy = 0;
        sumnumber = 0;
        for cir0 = 1:NumberOfBasement
            if cir0 ~= NumberOfBasement
                aoax = (b(cir0+1)-b(cir0))/(k(cir0)-k(cir0+1));
            else
                aoax = (b(cir0)-b(1))/(k(1)-k(cir0));
            end
            aoay = k(cir0).*aoax + b(cir0);
            
            sumx = sumx+aoax;
            sumy = sumy+aoay;
        end
        scatter(sumx/NumberOfBasement,sumy/NumberOfBasement,10,'b');
        data_aoa = [data_aoa;sumx/NumberOfBasement,sumy/NumberOfBasement];
     
    end
    csvwrite('AOA_PROCESSED.csv',data_aoa);
end
disp('Complete');

if TDOAEnable == 1 %TDOA SETUP--------------------------------------------------------------------------
disp('Printing TDOA Data');   
    c = 1;
    data = csvread('TDOA_RES.csv');
    lengthx = size(data);
    data_tdoa = [];
    base = [basex,basey];
    for i = 1:lengthx
        tdoa = data(i, :); 
        r = tdoa * c; 

        %LSM
        data_tdoa(i, :) = solve_tdoa(base, r);
    end
    scatter(data_tdoa(:,1), data_tdoa(:,2), 10, 'm');
end
csvwrite('TDOA_PROCESSED.csv',data_tdoa);
pause(1);
disp('Complete');


disp('Processing Data');
%--- Model 3 --%
figure(3);
grid on;
hold on;
xlabel('X/cm');
ylabel('Y/cm');
ylim([0,2000]);
axis equal;
set(gca,'ylim',[0 2000]);
set(gca,'xlim',[-1000 1000]);
% data_toa = %csvread('TOA_RES.csv');
% data_aoa = %csvread('AOA_RES.csv');
text1 = text(basex(1)+20, basey(1)+20, 'Base1', 'FontSize', 8, 'Color', 'red');
text2 = text(basex(2)+20, basey(2)+20, 'Base2', 'FontSize', 8, 'Color', 'red');
text3 = text(basex(3)+20, basey(3)+20, 'Base3', 'FontSize', 8, 'Color', 'red');
scatter(basex,basey,60,'g');
%TOA
if(TOAEnable == 1)
    x = data_toa(:,1);
    y = data_toa(:,2);
    fy1 = polyfit(x,y,1);
    y1_lsm = zeros( size(x));
    y1_lsm = polyval(fy1,x);
    plot(x,y1_lsm,'b-');
end
%AOA
if(AOAEnable == 1)
    x = data_aoa(:,1);
    y = data_aoa(:,2);
    fy1 = polyfit(x,y,1);
    y1_lsm = zeros( size(x));
    y1_lsm = polyval(fy1,x);
    plot(x,y1_lsm,'g-');
end
%TDOA
if(TDOAEnable == 1)
    x = data_tdoa(:,1);
    y = data_tdoa(:,2);
    fy1 = polyfit(x,y,1);
    % t = 0:0.1:lengthx;
    y1_lsm = zeros( size(x));
    y1_lsm = polyval(fy1,x);
    plot(x,y1_lsm,'m-');
end
disp('Complete');

checknum = 0;

if(TOAEnable == 1)
    data1 = data_toa;
    checknum = checknum+1;
else
    data1 = 0;
end
if(AOAEnable == 1)
    data2 = data_aoa;
    checknum = checknum+1;
else
    data2 = 0;
end
if(TDOAEnable == 1)
    data3 = data_tdoa;
    checknum = checknum+1;
else
    data3 = 0;
end

if(checknum == 0)
    disp('No Commend');
    return;
end

data = (data1+data2+data3)/checknum;
measured_x = data(:, 1)';
measured_y = data(:, 2)';
N = length(measured_x);

figure(4);
%plot(measured_x, measured_y, 'r.');
hold on;
grid on;


dt = 1; 
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]; 
H = [1 0 0 0; 0 0 1 0]; 
Q = 0.01 * eye(4); 
R = 0.1 * eye(2);

x = [measured_x(1); 0; measured_y(1); 0]; 
P = eye(4); 

filtered_x = zeros(1, N);
filtered_y = zeros(1, N);


for k = 1:N
 
    x = A * x;
    P = A * P * A' + Q;
    
    
    z = [measured_x(k); measured_y(k)]; 
    K = P * H' / (H * P * H' + R); 
    x = x + K * (z - H * x);
    P = (eye(4) - K * H) * P;
    
   
    filtered_x(k) = x(1); 
    filtered_y(k) = x(3); 
end

   data_original = csvread('OriginalLocation.csv');
   
   
   plot(data_original(:,1)/100,data_original(:,2)/100,'r', 'LineWidth', 1);
   
   
   %scatter((data_toa(:,1)+data_aoa(:,1))/200,(data_toa(:,2)+data_aoa(:,2))/200,10,'filled', 'MarkerFaceColor', 'b');
   scatter(data_tdoa(:,1)/100,data_tdoa(:,2)/100,10,'filled', 'MarkerFaceColor', 'b');
   %scatter(data_aoa(:,1)/100,data_aoa(:,2)/100,10,'filled', 'MarkerFaceColor', 'b');
   %scatter(data_tdoa(:,1)/100,data_tdoa(:,2)/100,10,'filled', 'MarkerFaceColor', 'b');
   
   plot(filtered_x/100, filtered_y/100, 'k', 'LineWidth', 1);
   
   legend('True trajectory', 'Least-square estimates', 'Tracking results');
   

% axis equal;
% set(gca,'ylim',[0 20]);
% set(gca,'xlim',[-10 10]);
% xlabel('X coordinate (m)', 'FontSize', 16);
% ylabel('Y coordinate (m)', 'FontSize', 16);
% disp('Complete');
% 
% 
% 
% figure(5);
% baseInfo = csvread('basementInfo.csv');
% NumberOfBasement = baseInfo(1,1);
% data_original = csvread('OriginalLocation.csv');
% data_processed1 = csvread('TDOA_PROCESSED.csv');
% %data_processed2 = csvread('AOA_PROCESSED.csv');
% %data_processed3 = csvread('TDOA_PROCESSED.csv');
% lengthx = size(data_original);
% 
% hold on;
% grid on;
% axis equal;
% set(gca,'ylim',[0 50], 'FontSize', 16);
% set(gca,'xlim',[0 100], 'FontSize', 16);
% datag = [];
% datab = [];
% cir = [];
% % 
% for cir0 = 1:lengthx
%     cir = [cir;cir0];
%     datag = [datag;sqrt((data_original(cir0,1)-filtered_x(cir0))^2+(data_original(cir0,2)-filtered_y(cir0))^2)];
%     datab = [datab;sqrt((data_original(cir0,1)-data_processed1(cir0,1))^2+(data_original(cir0,2)-data_processed1(cir0,2))^2)];
%     %datab = [datab;(sqrt((data_original(cir0,1)-data_processed1(cir0,1))^2+(data_original(cir0,2)-data_processed1(cir0,2))^2)+sqrt((data_original(cir0,1)-data_processed2(cir0,1))^2+(data_original(cir0,2)-data_processed2(cir0,2))^2))/2];
% end
% 
% scatter(cir(:,1),datab(:,1),10,'filled', 'MarkerFaceColor', 'b');
% plot(cir(:,1),datag(:,1),'r','LineWidth', 1);
% xlabel('Sample Point', 'FontSize', 16);
% ylabel('Location Error (cm)', 'FontSize', 16);
% legend('Least-square estimates', 'Tracking results')









function [distance] = DistantGet(xin1,yin1,xin2,yin2)
    distance = sqrt((yin2-yin1)^2+(xin2-xin1)^2);
end

function [value] = IntersectionCheck(xin1,yin1,r1,xin2,yin2,r2)
    distant = DistantGet(xin1,yin1,xin2,yin2);
    if distant > r1+r2
        value = 0;
    elseif distant == r1+r2
        value = 1;
    else
        value = 2;
    end
end

function [value] = unsign(vin)
    value = abs(vin);
end

function pos = solve_tdoa(BS, r)
    x0 = mean(BS, 1); 
    options = optimset('Display', 'off');
    pos = lsqnonlin(@(x) tdoa_residual(x, BS, r), x0, [], [], options);
end

function res = tdoa_residual(x, BS, r)
    d = vecnorm(BS - x, 2, 2); 
    d_diff = d(2:end) - d(1);
    res = d_diff - r(:);
end
