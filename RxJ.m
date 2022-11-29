

% love vars
a = -0.15;
b = 0.9;
y0 = [-1;1]; % initial conditions


% love formulas
% Q1
function dy = love_1(t,y,a,b)
dW = a*y(1) + b*y(2);
dE = 0;
dy = [dW;dE];
end
% Q2
function dy = love_2(t,y,a,b)
dW = a*y(1) + b*y(2);
dE = b*y(1) + a*y(2);
dy = [dW;dE];
end
% Q3
function dy = love_3(t,y,a,b)
dW = a*y(1) + b*y(2);
dE = -b*y(1) + -a*y(2);
dy = [dW;dE];
end
% Q4
function dy = love_4(t,y,a,b)
dW = a*y(1);
dE = b*y(1);
dy = [dW;dE];
end
% Q5
function dy = love_5(t,y,a,b)
dW = a*y(1) + b*y(2);
dE = b*y(1) + a*y(2);
dy = [dW;dE];
end

plove_1 = @(t,y) [a*y(1) + b*y(2); 0];
plove_2 = @(t,y) [a*y(1) + b*y(2); b*y(1) + a*y(2)];
plove_3 = @(t,y) [a*y(1) + b*y(2); -b*y(1) + -a*y(2)];
plove_4 = @(t,y) [a*y(1); b*y(1)];
%plove_5 = @(t,y) [a*y(1) + b*y(2); b*y(1) + a*y(2)];




%% doesn't change
% contruct graph
tstart = 0;
tfinal = 50;
tspan = [tstart tfinal];
% construct formulas
[t1,y1] = ode23(@(t,y) love_1(t,y,a,b),tspan,y0);
[t2,y2] = ode23(@(t,y) love_2(t,y,a,b),tspan,y0);
[t3,y3] = ode23(@(t,y) love_3(t,y,a,b),tspan,y0);
[t4,y4] = ode23(@(t,y) love_4(t,y,a,b),tspan,y0);
[t5,y5] = ode23(@(t,y) love_5(t,y,max(0,100-0.1*t),b),tspan,y0);

% display graph

figure
% 1
ax1 = subplot(3,2,1);
plot(ax1,t1,y1')
YLim = [-2 2];
XLabel.String = 'Time';
YLabel.String = 'Emotions';
Title.String = 'WallE & Eve''s';
%legend('WallE','Eve')
% 2
ax2 = subplot(3,2,2);
plot(ax2,t2,y2')
YLim = [-1.5 1.5];
XLabel.String = 'Time';
YLabel.String = 'Emotions';
Title.String = 'WallE & Eve''s relationship';
%legend('WallE','Eve')
% 3
ax3 = subplot(3,2,3);
plot(ax3,t3,y3')
YLim = [-1.5 1.5];
XLabel.String = 'Time';
YLabel.String = 'Emotions';
Title.String = 'WallE & Eve''s relationship';
%legend('WallE','Eve')
% 4
ax4 = subplot(3,2,4);
plot(ax4,t4,y4')
YLim = [-1.5 1.5];
XLabel.String = 'Time';
YLabel.String = 'Emotions';
Title.String = 'WallE & Eve''s relationship';
%legend('WallE','Eve')
% 5
ax5 = subplot(3,2,5);
plot(ax5,t5,y5')
YLim = [-1.5 1.5];
XLabel.String = 'Time';
YLabel.String = 'Emotions';
Title.String = 'WallE & Eve''s relationship';
%legend('WallE','Eve')


% portrait de phase

y1 = linspace(-10,10,20);
y2 = linspace(-10,10,20);
[uu,vv] = meshgrid(y2,y1);
init_time=0;
x11 = zeros(size(uu));
x12 = zeros(size(vv));
x21 = zeros(size(uu));
x22 = zeros(size(vv));
x31 = zeros(size(uu));
x32 = zeros(size(vv));
x41 = zeros(size(uu));
x42 = zeros(size(vv));
x51 = zeros(size(uu));
x52 = zeros(size(vv));

%a = -2e-1;
%b = 5e-1;
%love_1 = @(t,y) [a*y(2); b*y(1)];


% calculate data for portrait phases
for i = 1:numel(uu)
    %1
    Yder1 = plove_1(init_time,[uu(i); vv(i)]);
    x11(i) = Yder1(1);
    x12(i) = Yder1(2);
    %2
    Yder2 = plove_2(init_time,[uu(i); vv(i)]);
    x21(i) = Yder2(1);
    x22(i) = Yder2(2);
    %3
    Yder3 = plove_3(init_time,[uu(i); vv(i)]);
    x31(i) = Yder3(1);
    x32(i) = Yder3(2);
    %4
    Yder4 = plove_4(init_time,[uu(i); vv(i)]);
    x41(i) = Yder4(1);
    x42(i) = Yder4(2);
    %5
    %Yder5 = plove_5(init_time,[uu(i); vv(i)]);
    %x51(i) = Yder5(1);
    %x52(i) = Yder5(2);
end


% display portrait
figure 2
gca1 = subplot(3,2,1)
gca2 = subplot(3,2,2)
gca3 = subplot(3,2,3)
gca4 = subplot(3,2,4)
%gca5 = subplot(3,2,5)

%1
quiver(gca1,uu,vv,x11,x12,'r');
xlabel('Eve Emotions');
ylabel('WallE Emotions');
axis tight equal;
%2
quiver(gca2,uu,vv,x21,x22,'r');
xlabel('Eve Emotions');
ylabel('WallE Emotions');
axis tight equal;
%3
quiver(gca3,uu,vv,x31,x32,'r');
xlabel('Eve Emotions');
ylabel('WallE Emotions');
axis tight equal;
%4
quiver(gca4,uu,vv,x41,x42,'r');
xlabel('Eve Emotions');
ylabel('WallE Emotions');
axis tight equal;
%5
%quiver(gca5,uu,vv,x51,x52,'r');
%xlabel('Eve Emotions');
%ylabel('WallE Emotions');
%axis tight equal;





%y0_1 = [2;-1]; % initial conditions
%[t,y1] = ode23(@(t,y) plove_1(t,y),tspan,y0);
%figure(gcf)
%hold on
%plot(y1(:,2),y1(:,1),'b')








