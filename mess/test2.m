clear; clc; close all;


BS = [-1000, 200;
      -462, 926;
      0, 1800;
      345, 1882;
      1000, 1800;
      0, 200];

c = 1; 

a = 1.9;
b = 500; 

x_vals = linspace(-1000, 1000, 50); 
y_vals = a * x_vals + b; 
tdoa_data = zeros(length(x_vals), 5);

for i = 1:length(x_vals)
    target = [x_vals(i), y_vals(i)];
    
    % ???????????????
    dists = vecnorm(BS - target, 2, 2);
    
    r = dists(2:end) - dists(1);
    tdoa_data(i, :) = r / c;
end

csvwrite('TDOA_RES.csv', tdoa_data);

disp('TDOA data has been saved to TDOA_RES.csv');
