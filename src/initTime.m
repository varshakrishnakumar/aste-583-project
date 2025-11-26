function [ET0,LTM,LCM] = initTime()
% This function initializes the ephemieris time for the initial detection,
% the Lunar Targeting Maneuver, and the Lunar Correction Maneuver.

% Add Mice Path
addpath(genpath('mice\src\mice'));
addpath(genpath('mice\lib'));
addpath(genpath('mice\include'));
addpath(genpath('mice\generic_kernels'));

% Load DE440 Ephemeris
cspice_furnsh('mice\generic_kernels\de440s.bsp');

% Load Leapsecond Kernel
cspice_furnsh('mice\generic_kernels\naif0012.tls.pc.txt');

% Detection Date
date_detect = '1-DEC-2025 00:00:00.00 UTC'; % start date in UTC
ET0 = cspice_str2et(date_detect); % epoch time (seconds past J2000)

% LTM Date
date_LTM = '16-DEC-2025 00:00:00.00 UTC'; % start date in UTC
LTM = cspice_str2et(date_LTM); % epoch time (seconds past J2000)

% LCM Date
date_LCM = '25-DEC-2025 00:00:00.00 UTC'; % start date in UTC
LCM = cspice_str2et(date_LCM); % epoch time (seconds past J2000)

end