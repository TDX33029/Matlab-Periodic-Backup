%--   Version 1.0.2 DEBUG -- %
%------ Data Generate -------%

clear;
clc;
close all;
current_time = round(rem(now, 1) * 86400);
rng(current_time);

disp('Generate Application Setup');
%Config
SpreedSpeed = 3e8;
TrackingDelaySecond = 0.0;
isDisplayProcess = 0;
toa_error = 2;
aoa_error = 5;
tdoa_error = 2;

figure(1);
grid on;
hold on;
xlabel('Displacement X (m)');
ylabel('Displacement Y (m)');

baseInfo = csvread('basementInfo.csv');
NumberOfBasement = baseInfo(1,1);

axis equal;
set(gca,'ylim',[0 200]);
set(gca,'xlim',[-100 100]);
disp('Reading Info');
basex = csvread('basementInfo.csv',1,0,[1,0,NumberOfBasement,0]);
basey = csvread('basementInfo.csv',1,1,[1,1,NumberOfBasement,1]);
for cir0 = 1:NumberOfBasement
    label = sprintf('RU%d',cir0);
    text(basex(cir0)+2, basey(cir0)+2, label, 'FontSize', 8, 'Color', 'red');
end
scatter(basex,basey,30,'k');

%Get Tracking
h = drawfreehand('Color', 'blue', 'LineWidth', 2);

ObjectTracking = h.Position;

if size(ObjectTracking, 1) > 200
    step = ceil(size(ObjectTracking, 1) / 200);
    ObjectTracking = ObjectTracking(1:step:end, :);
end

num_points = 500;
if size(ObjectTracking, 1) > 2
    t = 1:size(ObjectTracking, 1);
    tt = linspace(1, size(ObjectTracking, 1), num_points);
    Objectx = spline(t, ObjectTracking(:,1), tt);
    Objecty = spline(t, ObjectTracking(:,2), tt);
    TrackingData = [Objectx', Objecty'];
else
    TrackingData = ObjectTracking;
end

delete(h);
for i = 1:size(TrackingData,1)-1
    plot(TrackingData(i:i+1,1), TrackingData(i:i+1,2), 'k', 'LineWidth', 1);
end
%Capture Done

%Save Original Data
csvwrite('Original_RES.csv',TrackingData);

scatter(TrackingData(1,1),TrackingData(1,2),15,'k');
textsetup = text(TrackingData(1,1)+5, TrackingData(1,2)+5, 'Starting Point', 'FontSize', 10, 'Color', 'red');
scatter(TrackingData(end,1),TrackingData(end,2),15,'k');
textsetup = text(TrackingData(end,1)+5, TrackingData(end,2)+5, 'Ending Point', 'FontSize', 10, 'Color', 'red');

%DataGenerate
r=2;

%TOA
disp('Generating TOA Data');
object = rectangle('Position',[TrackingData(1,1),TrackingData(1,2),2*r,2*r],'Curvature',[1,1],'FaceColor','b','EdgeColor','b');
textobject = text(TrackingData(1,1)+5,TrackingData(1,2)+5, 'Object', 'FontSize', 8, 'Color', 'b');
temp = [];
for cir1 = 1:num_points
    set(object, 'Position', [TrackingData(cir1,1)-r, TrackingData(cir1,2)-r, 2*r, 2*r]);
    set(textobject, 'Position', [TrackingData(cir1,1)-r+5, TrackingData(cir1,2)-r+5]);
    for cir0 = 1:NumberOfBasement
        distant_(cir0) = DistantGet(TrackingData(cir1,1),TrackingData(cir1,2),basex(cir0),basey(cir0));
    end
    
    %temp(end+1,1) = distant_(1).*normrnd(1,value1)/SpreedSpeed;
    temp(end+1,1) = (distant_(1)+normrnd(0,toa_error))/SpreedSpeed;
    
    for cir0 = 2:NumberOfBasement
       %temp(end,cir0) = distant_(cir0).*normrnd(1,value1)/SpreedSpeed;
       temp(end,cir0) = (distant_(cir0)+normrnd(0,toa_error))/SpreedSpeed;
    end
    
    if isDisplayProcess == 1
        pause(TrackingDelaySecond)
    end
    %drawnow;
end
% for cir0 = 1:NumberOfBasement
%     set(circle(cir0),'Position',[basex(cir0),basey(cir0), 0, 0]);
% end

csvwrite('TOA_RES.csv',temp);
disp('TOA Data Generated Completed');
pause(1);

%AOA
disp('Generating AOA Data');
object = rectangle('Position',[TrackingData(1,1),TrackingData(1,2),2*r,2*r],'Curvature',[1,1],'FaceColor','b','EdgeColor','b');
textobject = text(TrackingData(1,1)+5,TrackingData(1,2)+5, 'Object', 'FontSize', 8, 'Color', 'b');
temp = [];
for cir1 = 1:num_points
    set(object, 'Position', [TrackingData(cir1,1)-r, TrackingData(cir1,2)-r, 2*r, 2*r]);
    set(textobject, 'Position', [TrackingData(cir1,1)-r+5, TrackingData(cir1,2)-r+5]);
    %atan2(y1(spointx+cir1)-basey(cir0),spointx+cir1-basex(cir0));
    cir0 = 1;
    temp(end+1,1) = atan2(TrackingData(cir1,2)-basey(cir0),TrackingData(cir1,1)-basex(cir0))+normrnd(0,aoa_error*pi/180);
    for cir0 = 2:NumberOfBasement
        %temp(end,cir0) = atan2(y1(spointx+cir1)-basey(cir0),spointx+cir1-basex(cir0));
        temp(end,cir0) = atan2(TrackingData(cir1,2)-basey(cir0),TrackingData(cir1,1)-basex(cir0))+normrnd(0,aoa_error*pi/180);
    end

    if isDisplayProcess == 1
        pause(TrackingDelaySecond)
    end
    %drawnow;
end

csvwrite('AOA_RES.csv',temp);
disp('AOA Data Generated Completed');
pause(1);

%TDOA
TDOABasicRU = 1;    %参考基站:RU1
disp('Generating TDOA Data');
object = rectangle('Position',[TrackingData(1,1),TrackingData(1,2),2*r,2*r],'Curvature',[1,1],'FaceColor','b','EdgeColor','b');
textobject = text(TrackingData(1,1)+5,TrackingData(1,2)+5, 'Object', 'FontSize', 8, 'Color', 'b');
temp = [];
for cir1 = 1:num_points
    set(object, 'Position', [TrackingData(cir1,1)-r, TrackingData(cir1,2)-r, 2*r, 2*r]);
    set(textobject, 'Position', [TrackingData(cir1,1)-r+5, TrackingData(cir1,2)-r+5]);
    for cir0 = 2:NumberOfBasement
        length1 = DistantGet(TrackingData(cir1,1),TrackingData(cir1,2),basex(TDOABasicRU),basey(TDOABasicRU));
        length2 = DistantGet(TrackingData(cir1,1),TrackingData(cir1,2),basex(cir0),basey(cir0));
        % length0 = (length2-length1).*normrnd(1,value3)/SpreedSpeed;
        length0 = ((length2-length1)+normrnd(0,tdoa_error))/SpreedSpeed;
        if(cir0 == 2)
            temp(end+1,1) = length0;
        else
            temp(end,cir0-1) = length0;
        end
    end

    if isDisplayProcess == 1
        pause(TrackingDelaySecond)
    end
    %drawnow;
end
csvwrite('TDOA_RES.csv',temp);
disp('TDOA Data Generated Completed');

%Done
disp('Datas have been generated.');

%Function Unit
function [distance] = DistantGet(xin1,yin1,xin2,yin2)
    distance = sqrt((yin2-yin1)^2+(xin2-xin1)^2);
end
