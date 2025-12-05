function [] = plotGroundTrack(T,XrelEarth,params,st)
%This function plots the ground track of an earth-centered orbiter.
%   Inputs:
%       T = Time Vector
%       XrelEarth = Augmented State Vector Relative to Earth
%       params = project parameters
%       sc = project spacecraft parameters

% Constants
we = params.we; % Earth Rotation Rate in rad/s
GSlat = st.lat*180/pi; %station latitudes in deg
GSlong = st.long*180/pi; %station longitudes in deg
GST0 = params.GST0; %Initial Greenwich Meridian Time
% Position Vectors
rvec = XrelEarth(:,1:3);

% Determine Lat Long
lon = zeros(1,length(T));
lat = zeros(1,length(T));
for i = 1:1:length(T)
    % Rotation Matrix (EME to ECEF)
    R = [cos(we*T(i)+GST0) sin(we*T(i)+GST0) 0; -sin(we*T(i)+GST0) cos(we*T(i)+GST0) 0; 0 0 1];

    %Translate Vector
    rtemp = R*rvec(i,1:3)';
    
    lon(i) = atan2(rtemp(2),rtemp(1))*180/pi;
    lat(i) = asind(rtemp(3)/sqrt(rtemp(1)^2 + rtemp(2)^2 + rtemp(3)^2));
end

figure();
geoplot(lat,lon,'.')
hold on
geoplot(GSlat(1),GSlong(1),'.','MarkerSize',15,'Color','b')
geoplot(GSlat(2),GSlong(2),'.','MarkerSize',15,'Color','r')
geoplot(GSlat(3),GSlong(3),'.','MarkerSize',15,'Color','g')
geoplot(GSlat(4),GSlong(4),'.','MarkerSize',15,'Color','k')
title('Satellite Ground Track'); 
legend("Ground Track","Station 1","Station 2", "Station 3","Station 4")

end