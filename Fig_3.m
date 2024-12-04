clear;close all;clc;

minval = -1;
maxval = 1;
spaceval =  0.2;

[x, y] = meshgrid(minval:spaceval:maxval, minval:spaceval:maxval);

u_curl_free = x;  
v_curl_free = y; 

figure;
set(gcf,'position',[746 1 277 536])
subplot(2,1,1)
quiver(x, y, u_curl_free, v_curl_free, 'k','LineWidth',2);
xlabel('x');
ylabel('y');
axis off;
grid off;
drawnow

u_div_free = -y;  
v_div_free = x;  

subplot(2,1,2)
quiver(x, y, u_div_free, v_div_free, 'k','LineWidth',2);
xlabel('x');
ylabel('y');
axis off;
grid off;
drawnow