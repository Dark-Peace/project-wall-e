

% love formula
a = 2e-1;
b = 5e-1;
loveRJ = @(t,y) [-a*y(2); b*y(1)];



y0 = [-1;2]; % initial conditions


%% doesn't change
% contruct graph
tstart = 0;
tfinal = 50;
tspan = [tstart tfinal];
[t,y] = ode23(@(t,y) loveRJ(t,y),tspan,y0);
% display graph
plot(t,y')
ax = gca;
XLim = [0 50];
YLim = [-3 3];
XLabel.String = 'Time';
YLabel.String = 'Emotions';
String = 'Romeo & Juliet''s relationship';
legend('Romeo','Juliet')
