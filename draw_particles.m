path = '.';
files = dir([path '/*.csv']);  % automatically sorted by timestep

Lx=18e-6;
Ly=11e-6;
Lz=11e-6;
RADIUS=5e-6;

figure('pos',[100 100 600 600],'visible','off')
for j=1:length(files)
data = readtable([path '/' files(j).name]); % returns a table
R = mean(data.RADIUS);	% Particle radius

subplot(3,1,1);
t = 0:pi/20:2*pi;
yc = Ly/2;
zc = Lz/2;
y = yc + RADIUS*cos(t);
z = zc + RADIUS*sin(t);
plot(y,z,'color', 'k')
hold on
axis equal
title(['End view: ' num2str(j)])
xlabel('y'); ylabel('z')

subplot(3,1,2); % x,y plot
plot([0 Lx Lx 0],[Ly/2-RADIUS Ly/2-RADIUS Ly/2+RADIUS Ly/2+RADIUS],'color', 'k')  % outer square
hold on
axis equal
title(['Top view: ' num2str(j)])
xlabel('x'); ylabel('y')

subplot(3,1,3); % x,z plot
plot([0 Lx Lx 0],[Lz/2-RADIUS Lz/2-RADIUS Lz/2+RADIUS Lz/2+RADIUS],'color', 'k')  % outer square
hold on
axis equal
title(['Side view: ' num2str(j)])
xlabel('x'); ylabel('z')

for i=1:length(data.X)

subplot(3,1,1)
y = data.Y(i)+R*cos(t);
z = data.Z(i)+R*sin(t);
plot(y,z,'color','r')

subplot(3,1,2)
x = data.X(i)+R*cos(t);
y = data.Y(i)+R*sin(t);
plot(x,y,'color','r')
%xlim([0 Lx]);

subplot(3,1,3)
x = data.X(i)+R*cos(t);
z = data.Z(i)+R*sin(t);
plot(x,z,'color','r')
%xlim([0 Lx]);

end


subplot(3,1,1)
hold off

subplot(3,1,2)
hold off

subplot(3,1,3)
hold off
drawnow
saveas(gcf,['particles_' num2str(j) '.jpg']);
fprintf('saved plot for timestep %d\n',j);
end
