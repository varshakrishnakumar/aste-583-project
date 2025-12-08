function [observations] = computeObservations(T,X,X0fit,params,st) 
    w_e = params.we;
    GST0 = params.GST0; % Initial GST
    minEl = 10*(pi/180); %rad
    Re = params.Re; %km
    % Define Ground Stations
    GSlat = st.lat; %rad
    GSlat(4) = X0fit(8); %Update to OD solution for lat
    GSlong = st.long; %rad
    GSlat(4) = X0fit(9); %Update to OD solution for long

    measurements = [];
    for i = 1:1:length(T)
        % Step 1: Rotate r and range rate to topocentric frame
            LST1 = w_e*T(i)+GST0+GSlong(1); %LST at Station 1
            LST2 = w_e*T(i)+GST0+GSlong(2); %LST at Station 2
            LST3 = w_e*T(i)+GST0+GSlong(3); %LST at Station 3
            LST4 = w_e*T(i)+GST0+GSlong(4); %LST at Station 3
    
            %ECI to Topocentric Rotation Matrix for Station 1
            R1 = [sin(GSlat(1))*cos(LST1) sin(GSlat(1))*sin(LST1) -cos(GSlat(1));...
                  -sin(LST1) cos(LST1) 0;...
                  cos(GSlat(1))*cos(LST1) cos(GSlat(1))*sin(LST1) sin(GSlat(1))];
    
            %ECI to Topocentric Rotation Matrix for Station 2
            R2 = [sin(GSlat(2))*cos(LST2) sin(GSlat(2))*sin(LST2) -cos(GSlat(2));...
                  -sin(LST2) cos(LST2) 0;...
                  cos(GSlat(2))*cos(LST2) cos(GSlat(2))*sin(LST2) sin(GSlat(2))];
    
            %ECI to Topocentric Rotation Matrix for Station 3
            R3 = [sin(GSlat(3))*cos(LST3) sin(GSlat(3))*sin(LST3) -cos(GSlat(3));...
                  -sin(LST3) cos(LST3) 0;...
                  cos(GSlat(3))*cos(LST3) cos(GSlat(3))*sin(LST3) sin(GSlat(3))];

            %ECI to Topocentric Rotation Matrix for Station 4
            R4 = [sin(GSlat(4))*cos(LST4) sin(GSlat(4))*sin(LST4) -cos(GSlat(4));...
                  -sin(LST4) cos(LST4) 0;...
                  cos(GSlat(4))*cos(LST4) cos(GSlat(4))*sin(LST4) sin(GSlat(4))];
            
            %Rotate r and range rate
            rh1 = R1*X(i,1:3)';
            rh2 = R2*X(i,1:3)';
            rh3 = R3*X(i,1:3)';
            rh4 = R4*X(i,1:3)';
    
        % Step 2: Determine Range for each station
            range1 = rh1 - [0;0;Re];
            range2 = rh2 - [0;0;Re];
            range3 = rh3 - [0;0;Re];
            range4 = rh4 - [0;0;Re];
    
        % Step 3: Determine Range Rate for Each Station
            rangeRateGS1 = dot((X(i,1:3)-Re*[cos(LST1)*cos(GSlat(1)), sin(LST1)*cos(GSlat(1)), sin(GSlat(1))]),...
                            (X(i,4:6)-w_e*Re*[-sin(LST1)*cos(GSlat(1)), cos(LST1)*cos(GSlat(1)),0]))/norm(range1);
            rangeRateGS2 = dot((X(i,1:3)-Re*[cos(LST2)*cos(GSlat(2)), sin(LST2)*cos(GSlat(2)), sin(GSlat(2))]),...
                            (X(i,4:6)-w_e*Re*[-sin(LST2)*cos(GSlat(2)), cos(LST2)*cos(GSlat(2)),0]))/norm(range2);
            rangeRateGS3 = dot((X(i,1:3)-Re*[cos(LST3)*cos(GSlat(3)), sin(LST3)*cos(GSlat(3)), sin(GSlat(3))]),...
                            (X(i,4:6)-w_e*Re*[-sin(LST3)*cos(GSlat(3)), cos(LST3)*cos(GSlat(3)),0]))/norm(range3);
            rangeRateGS4 = dot((X(i,1:3)-Re*[cos(LST4)*cos(GSlat(4)), sin(LST4)*cos(GSlat(4)), sin(GSlat(4))]),...
                            (X(i,4:6)-w_e*Re*[-sin(LST4)*cos(GSlat(4)), cos(LST4)*cos(GSlat(4)),0]))/norm(range4);
    
        % Step 3: Determine Station Using El
            if(asin(range1(3)/norm(range1)) > minEl) %in view of station 1
                measurements = [measurements, [T(i);1;norm(range1);rangeRateGS1]];
            elseif(asin(range2(3)/norm(range2)) > minEl) %in view of station 2
                measurements = [measurements, [T(i);2;norm(range2);rangeRateGS2]];
            elseif(asin(range3(3)/norm(range3)) > minEl) %in view of station 3
                measurements = [measurements, [T(i);3;norm(range3);rangeRateGS3]];
            elseif(asin(range4(3)/norm(range4)) > minEl) %in view of station 3
                measurements = [measurements, [T(i);4;norm(range4);rangeRateGS4]];
            end
    
    end
    observations = measurements';
end