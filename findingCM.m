clear all;
close all;
clc;

%% Read in values, and divide into bins
positions = load('../Project3-build-desktop-Qt_4_6_2_in_PATH__System__Release/Random250Particlesxyz.txt');
x = positions(:,1);
y = positions(:,2);
z = positions(:,3);
vecLength = length(x);

[nElementsX, centersX] = hist(x,100);
[nElementsY, centersY] = hist(y,100);
[nElementsZ, centersZ] = hist(z,100);

%% Plot those histagrams, find those CM values
figure
hold on
hist(x,100);
fx = fit(centersX' ,nElementsX', 'gauss1')
plot(fx, centersX, nElementsX, 'g.')

figure
hold on
hist(y,100);
fy = fit(centersY' ,nElementsY', 'gauss1')
plot(fy, centersY, nElementsY, 'g.')

figure
hold on
hist(z,100);
fz = fit(centersZ' ,nElementsZ', 'gauss1')
plot(fz, centersZ, nElementsZ, 'g.')

% Lol, nope
r = zeros(1,vecLength);
for i=1:vecLength
    r(i) = sqrt( x(i)^2 + y(i)^2 + z(i)^2 );
end
[nElementsR, centersR] = hist(r,100);
figure
hold on
hist(r,100);
fr = fit(centersR', nElementsR', 'gauss1')
plot(fr, centersR, nElementsR)
%

%% Find density profile
% CM read from terminal
xCM = -3.216;
yCM = -3.219;
zCM =  3.494;

distance = zeros(1,vecLength);
for i=1:vecLength
    distance(i) = sqrt( (x(i)-xCM)^2 + (y(i)-yCM)^2 + (z(i)-zCM)^2 );
end
figure;
[n, bc] = hist(distance,100);
density = n./bc.^2;
plot(bc,density,'r.-')

n0 = 35;
r0 = 1.5;
r  = 0:0.1:bc(length(bc));
n = n0./(1+(r/r0).^4);
hold on 
plot(r,n,'k')

rho0 = 0.15;
rk = 15;
rho = rho0./( (r/rk).*(1+ (r/rk).^2) );
plot(r,rho)

rhofit = fit(centersR', nElementsR', 'rho0./( (r/rk).*(1+ (r/rk).^2) )')
plot

