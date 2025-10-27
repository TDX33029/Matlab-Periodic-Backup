%--  Version 1.0.9 DEBUG -- %

clear;
clc;
close all;

disp('Generate Application Setup');

value1 = 0.005;
value2 = 0.005;
value3 = 0.005;

figure(1);
grid on;
hold on;
final = [];
distant_ = [];
circle = [];
textdetail = [];

baseInfo = csvread('basementInfo.csv');
NumberOfBasement = baseInfo(1,1);

set(figure(1),'name','TOA');
textdetail0 = text(800,2000, 'Base Time: T','FontSize', 8, 'Color', 'red');
for cir0 = 1:NumberOfBasement
    textdetail(cir0) = text(800, 2000-cir0*50, '#UNDEFINE', 'FontSize', 8, 'Color', 'red');
end
% textdetail1 = text(800, 1950, ' ', 'FontSize', 8, 'Color', 'red');
% textdetail2 = text(800, 1900, ' ', 'FontSize', 8, 'Color', 'red');
% textdetail3 = text(800, 1850, ' ', 'FontSize', 8, 'Color', 'red');
% textdetail4 = text(800, 1850, ' ', 'FontSize', 8, 'Color', 'red');
% axis([0,2000,0,2000]);
% xlim([0,2000]);
% ylim([0,2000]);
axis equal;
set(gca,'ylim',[0 2000]);
set(gca,'xlim',[-1000 1000]);
disp('Reading Info');
basex = csvread('basementInfo.csv',1,0,[1,0,NumberOfBasement,0]);
basey = csvread('basementInfo.csv',1,1,[1,1,NumberOfBasement,1]);
for cir0 = 1:NumberOfBasement
    label = sprintf('Base%d',cir0);
    text(basex(cir0)+20, basey(cir0)+20, label, 'FontSize', 8, 'Color', 'red');
end
disp('Capturing Location');
scatter(basex,basey,60,'g');
[spointx,spointy] = ginput(1);
scatter(spointx,spointy,60,'r','fill');
textsetup = text(spointx+20, spointy+20, 'Setup', 'FontSize', 8, 'Color', 'red');
[epointx,epointy] = ginput(1);
scatter(epointx,epointy,60,'r');
textend = text(epointx+20, epointy+20, 'End', 'FontSize', 8, 'Color', 'red');
plot([spointx,epointx],[spointy,epointy],'r');
k1 = (epointy-spointy)/(epointx-spointx);
y1 = @(x)k1.*x+spointy-k1.*spointx;
r = 10;
for cir0 = 1:NumberOfBasement
    distant_(cir0) = DistantGet(spointx,spointy,basex(cir0),basey(cir0));
    circle(cir0) = rectangle('Position', [basex(cir0)-distant_(cir0),basey(cir0)-distant_(cir0), 2*distant_(cir0), 2*distant_(cir0)],...
          'Curvature', [1, 1],...
          'EdgeColor', 'b',...
          'LineWidth', 0.2);
end
% distant1 = DistantGet(spointx,spointy,basex(1),basey(1));
% distant2 = DistantGet(spointx,spointy,basex(2),basey(2));
% distant3 = DistantGet(spointx,spointy,basex(3),basey(3));

% circle1 = rectangle('Position', [basex(1)-distant1,basey(1)-distant1, 2*distant1, 2*distant1],...
%           'Curvature', [1, 1],...
%           'EdgeColor', 'b',...
%           'LineWidth', 0.2);
% circle2 = rectangle('Position', [basex(2)-distant2,basey(2)-distant2, 2*distant2, 2*distant2],...
%           'Curvature', [1, 1],...
%           'EdgeColor', 'b',...
%           'LineWidth', 0.2);
% circle3 = rectangle('Position', [basex(3)-distant3,basey(3)-distant3, 2*distant3, 2*distant3],...
%           'Curvature', [1, 1],...
%           'EdgeColor', 'b',...
%           'LineWidth', 0.2);
object = rectangle('Position',[spointx-r,spointy-r,2*r,2*r],'Curvature',[1,1],'FaceColor','b','EdgeColor','b');
textobject = text(spointx+20, spointy+20, 'Object', 'FontSize', 8, 'Color', 'b');

if epointx-spointx>=0
    lengthx = epointx-spointx;
    Forward = 0;
else
    lengthx = spointx-epointx;
    Forward = 1;
end

disp('Generating TOA Data');
%--- TOA ---%
temp = [];
for cir1 = 1:lengthx
    if Forward == 0
        set(object, 'Position', [spointx+cir1-r, y1(spointx+cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx+cir1-r+20, y1(spointx+cir1)-r+20]);
        for cir0 = 1:NumberOfBasement
            distant_(cir0) = DistantGet(spointx+cir1,y1(spointx+cir1),basex(cir0),basey(cir0));
        end
%         distant1 = DistantGet(spointx+cir1,y1(spointx+cir1),basex(1),basey(1));
%         distant2 = DistantGet(spointx+cir1,y1(spointx+cir1),basex(2),basey(2));
%         distant3 = DistantGet(spointx+cir1,y1(spointx+cir1),basex(3),basey(3));
    elseif Forward == 1
        set(object, 'Position', [spointx-cir1-r, y1(spointx-cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx-cir1-r+20, y1(spointx-cir1)-r+20]);
        for cir0 = 1:NumberOfBasement
            distant_(cir0) = DistantGet(spointx-cir1,y1(spointx-cir1),basex(cir0),basey(cir0));
        end
%         distant1 = DistantGet(spointx-cir1,y1(spointx-cir1),basex(1),basey(1));
%         distant2 = DistantGet(spointx-cir1,y1(spointx-cir1),basex(2),basey(2));
%         distant3 = DistantGet(spointx-cir1,y1(spointx-cir1),basex(3),basey(3));
    else
        disp('Error');
        return   %Error Exit
    end
    
    
    %set(circle(1),'Position',[basex(1)-distant_(1),basey(1)-distant_(1), 2*distant_(1), 2*distant_(1)]); 
    temp(end+1,1) = distant_(1).*normrnd(1,value1);
    
    for cir0 = 2:NumberOfBasement
       %set(circle(cir0),'Position',[basex(cir0)-distant_(cir0),basey(cir0)-distant_(cir0), 2*distant_(cir0), 2*distant_(cir0)]); 
       textdetail_ = sprintf('Base1 Receive Time: %.2f',distant_(cir0));
       %set(textdetail(cir0), 'String', textdetail_);
       temp(end,cir0) = distant_(cir0).*normrnd(1,value1);
    end
%     set(circle1,'Position',[basex(1)-distant1,basey(1)-distant1, 2*distant1, 2*distant1]);
%     set(circle2,'Position',[basex(2)-distant2,basey(2)-distant2, 2*distant2, 2*distant2]);
%     set(circle3,'Position',[basex(3)-distant3,basey(3)-distant3, 2*distant3, 2*distant3]);
%     textdetail1_ = sprintf('Base1 Receive Time: %.2f',distant1);
%     textdetail2_ = sprintf('Base2 Receive Time: %.2f',distant2);
%     textdetail3_ = sprintf('Base3 Receive Time: %.2f',distant3);
    
%     set(textdetail1, 'String', textdetail1_);
%     set(textdetail2, 'String', textdetail2_);
%     set(textdetail3, 'String', textdetail3_);
    drawnow;
    % pause(0.001);
end
for cir0 = 1:NumberOfBasement
    set(circle(cir0),'Position',[basex(cir0),basey(cir0), 0, 0]);
end
% set(circle1,'Position',[basex(1),basey(1), 0, 0]);
% set(circle2,'Position',[basex(2),basey(2), 0, 0]);
% set(circle3,'Position',[basex(3),basey(3), 0, 0]);
csvwrite('TOA_RES.csv',temp);
disp('Complete');
pause(1);

disp('Generating AOA Data');
%--- AOA ---%
k = [];
temp = [];
handler = [];

object = rectangle('Position',[spointx-r,spointy-r,2*r,2*r],'Curvature',[1,1],'FaceColor','b','EdgeColor','b');
textobject = text(spointx+20, spointy+20, 'Object', 'FontSize', 8, 'Color', 'b');
for cir1 = 1:lengthx
    if Forward == 0
        set(object, 'Position', [spointx+cir1-r, y1(spointx+cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx+cir1-r+20, y1(spointx+cir1)-r+20]);
        k(1) = (y1(spointx+cir1)-basey(1))/(spointx+cir1-basex(1));%tan(temp(end+1,1));
        temp(end+1,1) = atan(k(1)).*normrnd(1,value2);
        for cir0 = 2:NumberOfBasement
           k(cir0) =  (y1(spointx+cir1)-basey(cir0))/(spointx+cir1-basex(cir0));%tan(temp(end,cir0));
           temp(end,cir0) =  atan(k(cir0)).*normrnd(1,value2);
        end
    elseif Forward == 1
        set(object, 'Position', [spointx-cir1-r, y1(spointx-cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx-cir1-r+20, y1(spointx-cir1)-r+20]);
        k(1) = (y1(spointx-cir1)-basey(1))/(spointx-cir1-basex(1));%tan(temp(end,1));
        temp(end+1,1) = atan(k(1)).*normrnd(1,value2);
        for cir0 = 2:NumberOfBasement
           k(cir0) = (y1(spointx-cir1)-basey(cir0))/(spointx-cir1-basex(cir0));%tan(temp(end,cir0));
           temp(end,cir0) =  atan(k(cir0)).*normrnd(1,value2);
        end
    else
        disp('Error');
        return
    end
    for cir0 = 1:NumberOfBasement
        
        y2 = @(x)k(cir0).*(x-basex(cir0))+basey(cir0);
        %handler(cir0) = fplot(y2, [-1000, 1000], 'b', 'LineWidth', 0.01);   
    end
    pause(0.0);
    %delete(handler);
end
csvwrite('AOA_RES.csv',temp);
disp('Complete');
pause(1);

disp('Generating TDOA Data');   %-----------TDOA-------------------------------------------------
temp = [];
for cir1  = 1:lengthx
    if Forward == 0
        set(object, 'Position', [spointx+cir1-r, y1(spointx+cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx+cir1-r+20, y1(spointx+cir1)-r+20]);
        for cir0 = 2:NumberOfBasement
            length1 = DistantGet(spointx+cir1,y1(spointx+cir1),basex(1),basey(1));
            length2 = DistantGet(spointx+cir1,y1(spointx+cir1),basex(cir0),basey(cir0));
            %scatter(spointx+cir1,y1(spointx+cir1),30,'b');
            length0 = (length2-length1).*normrnd(1,value3);
            if(cir0 == 2)
                temp(end+1,1) = length0;
            else
                temp(end,cir0-1) = length0;
            end
        end
        
    elseif Forward == 1
        set(object, 'Position', [spointx-cir1-r, y1(spointx-cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx-cir1-r+20, y1(spointx-cir1)-r+20]);
        for cir0 = 2:NumberOfBasement
            length1 = DistantGet(spointx-cir1,y1(spointx-cir1),basex(1),basey(1));
            length2 = DistantGet(spointx-cir1,y1(spointx-cir1),basex(cir0),basey(cir0));
            %scatter(spointx-cir1,y1(spointx-cir1),30,'b');
            length0 = (length2-length1).*normrnd(1,value3);
            if(cir0 == 2)
                temp(end+1,1) = length0;
            else
                temp(end,cir0-1) = length0;
            end
        end
        
    else
        disp('Error');
        return
    end
    pause(0.0);
end
csvwrite('TDOA_RES.csv',temp);
disp('Complete');


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
%

list = [];
baseInfo = csvread('basementInfo.csv');
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
% axis([0,2000,0,2000]);
% xlim([0,2000]);
ylim([0,2000]);
axis equal;
set(gca,'ylim',[0 2000]);
set(gca,'xlim',[-1000 1000]);


if TOAEnable == 1   %TOA SETUP---------------------------------------------------------------------
disp('Printing TOA Data');
    data = csvread('TOA_RES.csv');
    for cir0 = 1:NumberOfBasement
        label = sprintf('Base%d',cir0);
        text(basex(cir0)+20, basey(cir0)+20, label, 'FontSize', 8, 'Color', 'red');
    end
    % text1 = text(basex(1)+20, basey(1)+20, 'Base1', 'FontSize', 8, 'Color', 'red');
    % text2 = text(basex(2)+20, basey(2)+20, 'Base2', 'FontSize', 8, 'Color', 'red');
    % text3 = text(basex(3)+20, basey(3)+20, 'Base3', 'FontSize', 8, 'Color', 'red');
    scatter(spointx,spointy,60,'r','fill');
    textsetup = text(spointx+20, spointy+20, 'Setup', 'FontSize', 8, 'Color', 'red');
    scatter(epointx,epointy,60,'r');
    textend = text(epointx+20, epointy+20, 'End', 'FontSize', 8, 'Color', 'red');
    plot([spointx,epointx],[spointy,epointy],'r');
    scatter(basex,basey,60,'g');

    for cir1 = 1:NumberOfBasement
        circle0 = rectangle('Position', [basex(cir1)-data(1,cir1),basey(cir1)-data(1,cir1),1,1],...
              'Curvature', [1, 1],...
              'EdgeColor', 'b',...
              'LineWidth', 0.2);
    end
    % circle1 = rectangle('Position', [basex(1)-data(1,1),basey(1)-data(1,1),1,1],...
    %           'Curvature', [1, 1],...
    %           'EdgeColor', 'b',...
    %           'LineWidth', 0.2);
    % circle2 = rectangle('Position', [basex(2)-data(1,2),basey(2)-data(1,2),1,1],...
    %           'Curvature', [1, 1],...
    %           'EdgeColor', 'b',...
    %           'LineWidth', 0.2);
    % circle3 = rectangle('Position', [basex(3)-data(1,3),basey(3)-data(1,3),1,1],...
    %           'Curvature', [1, 1],...
    %           'EdgeColor', 'b',...
    %           'LineWidth', 0.2);

%     point12 = [0,0];
%     point13 = [0,0];
%     point23 = [0,0];

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
            sum(1) = sum(1) +  point(1);
            sum(2) = sum(2) +  point(2);
%             scatter(point(1),point(2),10,'b');     %Print TOA Details
        end
        scatter(sum(1)/NumberOfBasement,sum(2)/NumberOfBasement,10,'r');
        data_toa = [data_toa;sum(1)/NumberOfBasement,sum(2)/NumberOfBasement];
    %         
    %     temp1 = IntersectionCheck(basex(1),basey(1),data(cir1,1),basex(2),basey(2),data(cir1,2));
    %     temp2 = IntersectionCheck(basex(1),basey(1),data(cir1,1),basex(3),basey(3),data(cir1,3));
    %     temp3 = IntersectionCheck(basex(2),basey(2),data(cir1,2),basex(3),basey(3),data(cir1,3));
    %     1-2
    %     if( temp1 == 0 || temp1 == 1)
    %         point12(1) = basex(1)+((basex(2)-basex(1))/(data(cir1,1)+data(cir1,2)))*data(cir1,1);
    %         point12(2) = basey(1)+((basey(2)-basey(1))/(data(cir1,1)+data(cir1,2)))*data(cir1,1);
    %     elseif(temp1 == 2)
    %         d = DistantGet(basex(1),basey(1),basex(2),basey(2));
    %         m = (data(cir1,1)^2+d^2-data(cir1,2)^2)/(2*d);
    %         h = sqrt(data(cir1,1)^2-m^2);
    %         point12_m(1) = basex(1)+((basex(2)-basex(1))/d)*m;
    %         point12_m(2) = basey(1)+((basey(2)-basey(1))/d)*m;
    %         point12_t(1) = point12_m(1)+h*(basey(2)-basey(1))/d;
    %         point12_t(2) = point12_m(2)-h*(basex(2)-basex(1))/d;
    %         point12_t(3) = point12_m(1)-h*(basey(2)-basey(1))/d;
    %         point12_t(4) = point12_m(2)+h*(basex(2)-basex(1))/d;
    %         
    %         d3(1) =  DistantGet(point12_t(1),point12_t(2),basex(3),basey(3));
    %         d3(2) =  DistantGet(point12_t(3),point12_t(4),basex(3),basey(3));
    %         if(unsign(data(cir1,3) - d3(1)) > unsign(data(cir1,3) - d3(2)))
    %             point12(1) = point12_t(3);
    %             point12(2) = point12_t(4);
    %         else
    %             point12(1) = point12_t(1);
    %             point12(2) = point12_t(2);
    %         end
    %     else
    %         disp('Exit12');
    %         return %exit
    %     end
    %     1-3
    %     if( temp2 == 0 || temp2 == 1)
    %         point13(1) = basex(1)+((basex(3)-basex(1))/(data(cir1,1)+data(cir1,3)))*data(cir1,1);
    %         point13(2) = basey(1)+((basey(3)-basey(1))/(data(cir1,1)+data(cir1,3)))*data(cir1,1);
    %     elseif(temp2 == 2)
    %         d = DistantGet(basex(1),basey(1),basex(3),basey(3));
    %         m = (data(cir1,1)^2+d^2-data(cir1,3)^2)/(2*d);
    %         h = sqrt(data(cir1,1)^2-m^2);
    %         point13_m(1) = basex(1)+((basex(3)-basex(1))/d)*m;
    %         point13_m(2) = basey(1)+((basey(3)-basey(1))/d)*m;
    %         point13_t(1) = point13_m(1)+h*(basey(3)-basey(1))/d;
    %         point13_t(2) = point13_m(2)-h*(basex(3)-basex(1))/d;
    %         point13_t(3) = point13_m(1)-h*(basey(3)-basey(1))/d;
    %         point13_t(4) = point13_m(2)+h*(basex(3)-basex(1))/d;
    %         
    %         d3(1) =  DistantGet(point13_t(1),point13_t(2),basex(2),basey(2));
    %         d3(2) =  DistantGet(point13_t(3),point13_t(4),basex(2),basey(2));
    %         if(unsign(data(cir1,2) - d3(1)) > unsign(data(cir1,2) - d3(2)))
    %             point13(1) = point13_t(3);
    %             point13(2) = point13_t(4);
    %         else
    %             point13(1) = point13_t(1);
    %             point13(2) = point13_t(2);
    %         end
    %     else
    %         disp('Exit13');
    %         return %exit
    %     end
    %     2-3
    %     if( temp3 == 0 || temp3 == 1)
    %         point23(1) = basex(2)+((basex(3)-basex(2))/(data(cir1,2)+data(cir1,3)))*data(cir1,2);
    %         point23(2) = basey(2)+((basey(3)-basey(2))/(data(cir1,2)+data(cir1,3)))*data(cir1,2);
    %     elseif(temp3 == 2)
    %         d = DistantGet(basex(2),basey(2),basex(3),basey(3));
    %         m = (data(cir1,2)^2+d^2-data(cir1,3)^2)/(2*d);
    %         h = sqrt(data(cir1,2)^2-m^2);
    %         point23_m(1) = basex(2)+((basex(3)-basex(2))/d)*m;
    %         point23_m(2) = basey(2)+((basey(3)-basey(2))/d)*m;
    %         point23_t(1) = point23_m(1)+h*(basey(3)-basey(2))/d;
    %         point23_t(2) = point23_m(2)-h*(basex(3)-basex(2))/d;
    %         point23_t(3) = point23_m(1)-h*(basey(3)-basey(2))/d;
    %         point23_t(4) = point23_m(2)+h*(basex(3)-basex(2))/d;
    %         
    %         d3(1) =  DistantGet(point23_t(1),point23_t(2),basex(1),basey(1));
    %         d3(2) =  DistantGet(point23_t(3),point23_t(4),basex(1),basey(1));
    %         if(unsign(data(cir1,1) - d3(1)) > unsign(data(cir1,1) - d3(2)))
    %             point23(1) = point23_t(3);
    %             point23(2) = point23_t(4);
    %         else
    %             point23(1) = point23_t(1);
    %             point23(2) = point23_t(2);
    %         end
    %     else
    %         disp('Exit23');
    %         return %exit
    %     end

    %     %--
    %     scatter(point12(1),point12(2),10,'b');
    %     scatter(point13(1),point13(2),10,'b');
    %     scatter(point23(1),point23(2),10,'b');
    %     scatter((point12(1)+point13(1)+point23(1))/3,(point12(2)+point13(2)+point23(2))/3,10,'r');
    %     
    %     set(circle1,'Position',[basex(1)-data(cir1,1),basey(1)-data(cir1,1), 2*data(cir1,1), 2*data(cir1,1)]);
    %     set(circle2,'Position',[basex(2)-data(cir1,2),basey(2)-data(cir1,2), 2*data(cir1,2), 2*data(cir1,2)]);
    %     set(circle3,'Position',[basex(3)-data(cir1,3),basey(3)-data(cir1,3), 2*data(cir1,3), 2*data(cir1,3)]);
    %     final = [final;(point12(1)+point13(1)+point23(1))/3,(point12(2)+point13(2)+point23(2))/3];
        %pause(0.01);
    end
    for cir0 = 1:NumberOfBasement
        set(circle(cir0),'Position',[basex(cir0),basey(cir0), 0, 0]);
    end
    % set(circle1,'Position',[basex(1),basey(1), 0, 0]);
    % set(circle2,'Position',[basex(2),basey(2), 0, 0]);
    % set(circle3,'Position',[basex(3),basey(3), 0, 0]);
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
    for cir1 = 1:lengthx
        for cir0 = 1:NumberOfBasement
            k(cir0) = tan(data(cir1,cir0));
            b(cir0) = basey(cir0) - k(cir0).*basex(cir0);
        end
        for cir0 = 1:NumberOfBasement
            y3 = @(x)k(cir0).*(x-basex(cir0))+basey(cir0);
            %handler(cir0) = fplot(y3, [-1000, 1000], 'b', 'LineWidth', 0.01);
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
        
        pause(0.0);
        %delete(handler);
    end
    
end
disp('Complete');

if TDOAEnable == 1 %TDOA SETUP--------------------------------------------------------------------------
disp('Printing TDOA Data');   
    c = 1;
    data = csvread('TDOA_RES.csv');
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

pause(1);
disp('Complete');


disp('Processing Data');
%--- Model 3 --%
figure(3);
grid on;
hold on;
% axis([0,2000,0,2000]);
% xlim([0,2000]);
ylim([0,2000]);
axis equal;
set(gca,'ylim',[0 2000]);
set(gca,'xlim',[-1000 1000]);
% data_toa = %csvread('TOA_RES.csv');
% data_aoa = %csvread('AOA_RES.csv');
text1 = text(basex(1)+20, basey(1)+20, 'Base1', 'FontSize', 8, 'Color', 'red');
text2 = text(basex(2)+20, basey(2)+20, 'Base2', 'FontSize', 8, 'Color', 'red');
text3 = text(basex(3)+20, basey(3)+20, 'Base3', 'FontSize', 8, 'Color', 'red');
scatter(spointx,spointy,60,'r','fill');
textsetup = text(spointx+20, spointy+20, 'Setup', 'FontSize', 8, 'Color', 'red');
scatter(epointx,epointy,60,'r');
textend = text(epointx+20, epointy+20, 'End', 'FontSize', 8, 'Color', 'red');
plot([spointx,epointx],[spointy,epointy],'r');
scatter(basex,basey,60,'g');
%TOA
if(TOAEnable == 1)
    x = data_toa(:,1);
    y = data_toa(:,2);
    fy1 = polyfit(x,y,1);
    % t = 0:0.1:lengthx;
    y1_lsm = zeros( size(x));
    y1_lsm = polyval(fy1,x);
    plot(x,y1_lsm,'b-');
end
%AOA
if(AOAEnable == 1)
    x = data_aoa(:,1);
    y = data_aoa(:,2);
    fy1 = polyfit(x,y,1);
    % t = 0:0.1:lengthx;
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


% Read data

% TOAEnable = 1;
% AOAEnable = 1;
% TDOAEnable = 1;

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


plot(measured_x, measured_y, 'r.');
hold on;
grid on;

% Kalman filter parameters
dt = 1; % Time step size (assuming data points are equally spaced)
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1]; % State transition matrix
H = [1 0 0 0; 0 0 1 0]; % Observation matrix
Q = 0.01 * eye(4); %Process noise covariance matrix
R = 0.1 * eye(2); % Observed noise covariance matrix

% Initialization state and covariance matrix
x = [measured_x(1); 0; measured_y(1); 0]; % Initial state [x, vx, y, vy]
P = eye(4); % Initial covariance matrix

% Stores the filtered state
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


plot(filtered_x, filtered_y, 'g', 'LineWidth', 2);

axis equal;
set(gca,'ylim',[0 2000]);
set(gca,'xlim',[-1000 1000]);


% % RTS????
% smoothed_x = zeros(1, N);
% smoothed_y = zeros(1, N);
% 
% % ?????????????
% x_smooth = x; % ?????????????????????
% P_smooth = P; % ???????????????????????????
% 
% % ???????????
% smoothed_x(N) = x_smooth(1);
% smoothed_y(N) = x_smooth(3);
% 
% % ??????
% for k = N-1:-1:1
%     % ?????????????????????
%     x_k = [filtered_x(k); 0; filtered_y(k); 0]; % ????
%     P_k = P; % ???????
%     
%     % ????
%     x_pred = A * x_smooth;
%     P_pred = A * P_smooth * A' + Q;
%     
%     % ????
%     C = P_k * A' / P_pred;
%     
%     % ??????
%     x_smooth = x_k + C * (x_smooth - x_pred);
%     P_smooth = P_k + C * (P_smooth - P_pred) * C';
%     
%     % ????????
%     smoothed_x(k) = x_smooth(1);
%     smoothed_y(k) = x_smooth(3);
% end
% 
% % ????????
% plot(smoothed_x, smoothed_y, 'm', 'LineWidth', 2);
% legend('????', '?????', '?????');
% title('??????????????');
% xlabel('X');
% ylabel('Y');
% grid on;
disp('Complete');

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

% TDOA ??????
function pos = solve_tdoa(BS, r)
    % ??????????
    x0 = mean(BS, 1); % ?????
    options = optimset('Display', 'off');
    pos = lsqnonlin(@(x) tdoa_residual(x, BS, r), x0, [], [], options);
end

function res = tdoa_residual(x, BS, r)
    % ?? TDOA ??
    d = vecnorm(BS - x, 2, 2); % ??????????????
    d_diff = d(2:end) - d(1);  % ?????1????
    res = d_diff - r(:);       % ????
end
