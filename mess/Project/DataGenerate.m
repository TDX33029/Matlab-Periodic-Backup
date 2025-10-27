%--  Version 1.0.10 DEBUG -- %
%------  Data Generate-------%

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

    elseif Forward == 1
        set(object, 'Position', [spointx-cir1-r, y1(spointx-cir1)-r, 2*r, 2*r]);
        set(textobject, 'Position', [spointx-cir1-r+20, y1(spointx-cir1)-r+20]);
        for cir0 = 1:NumberOfBasement
            distant_(cir0) = DistantGet(spointx-cir1,y1(spointx-cir1),basex(cir0),basey(cir0));
        end

    else
        disp('Error');
        return   %Error Exit
    end
    
   temp(end+1,1) = distant_(1).*normrnd(1,value1);
    
    for cir0 = 2:NumberOfBasement
       textdetail_ = sprintf('Base1 Receive Time: %.2f',distant_(cir0));
       temp(end,cir0) = distant_(cir0).*normrnd(1,value1);
    end

    drawnow;
end
for cir0 = 1:NumberOfBasement
    set(circle(cir0),'Position',[basex(cir0),basey(cir0), 0, 0]);
end

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
        k(1) = (y1(spointx-cir1)-basey(1))/(spointx-cir1-basex(1));
        temp(end+1,1) = atan(k(1)).*normrnd(1,value2);
        for cir0 = 2:NumberOfBasement
           k(cir0) = (y1(spointx-cir1)-basey(cir0))/(spointx-cir1-basex(cir0));
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
disp('Datas have been generated.');
%Function Unit

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


