clear
clc
[params,sc,st,X0,P0] = projectInit();
[ET0,LTM,LCM] = initTime();
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
step = 1800;
X0preLTM = [X0;reshape(eye(10),100,1)];
X0preLTM(7) = 0.6996824;
[TpreLTM,XpreLTM] = ode45(@projectDynamics,ET0:step:LTM,X0preLTM,options,params,sc);
X0postLTM = XpreLTM(end,:)';
deltaVLTM = [0.6931075; -0.8462091; 0.0956979];
X0postLTM(4:6) = X0postLTM(4:6) + deltaVLTM;
[TpostLTM,XpostLTM] = ode45(@projectDynamics,LTM:step:LCM,X0postLTM,options,params,sc);
X0postLCM = XpostLTM(end,:)'; %Don't apply to reference trajectory
TOF = 227*24*3600;
[TpostLCM,XpostLCM] = ode45(@projectDynamics,LCM:step:LCM+TOF,X0postLCM,options,params,sc);
Xtotal = [XpreLTM(:,:);XpostLTM(2:end,:);XpostLCM(2:end,:)];
Ttotal = [TpreLTM;TpostLTM(2:end);TpostLCM(2:end)];
deltaVLTM_nom = [0.6931075; -0.8462091; 0.0956979];
xEarth = cspice_spkezr('Earth', Ttotal.', 'ECLIPJ2000','NONE','Sun');  % 6xN
xEarth = xEarth.';   % Nx6
xMoon  = cspice_spkezr('Moon',  Ttotal.', 'ECLIPJ2000','NONE','Earth'); % 6xN
xMoon  = xMoon.';    % Nx6
XrelEarth = Xtotal(:,1:6) - xEarth;
opts = struct();
opts.maxIter    = 11;
opts.earlyMode  = "drop";
opts.earlyDays  = 2;
opts.earlyScale = 10;
out = runBatchLSQ('0to14','ekf_k',opts);
deltaX0 = out.X_est - X0;
P_final = out.P_est;
disp('deltaX0 ='); disp(deltaX0.');
disp('P_final ='); disp(P_final);
X0fit = X0 + deltaX0;
X0fit(7) = 0.6996824;
Pfit = P_final;
ntrials = 1000;
LCM_th = 30*params.tday;
indexLCM    = floor((LCM-ET0)/step) + 1;
indexLCM_th = indexLCM + LCM_th/step;
STM_LCM_Target = reshape(Xtotal(indexLCM_th,11:end),10,10) * ...
                 inv(reshape(Xtotal(indexLCM,11:end),10,10));
deltaVLTMmag = zeros(ntrials,1);
deltaVLCMmag = zeros(ntrials,1);
posMagError  = zeros(ntrials,1);
velMagError  = zeros(ntrials,1);
Rerror       = zeros(ntrials,1);
Verror       = zeros(ntrials,1);
for i = 1:ntrials
    stateDeviation      = randn(10,1).*sqrt(diag(Pfit));
    stateDeviation(7)   = 0;
    X0samp              = X0fit + stateDeviation;
    X0sampPreLTM = X0samp;
    [TpreLTM,XpreLTM] = ode45(@projectDynamics, ET0:step:LTM, ...
                              X0sampPreLTM, options, params, sc);
    deltaVLTM = [0.6931075; -0.8462091; 0.0956979] + randn(3,1)*0.005;
    deltaVLTMmag(i) = norm(deltaVLTM);
    X0sampPostLTM      = XpreLTM(end,:)';   % still 10x1
    X0sampPostLTM(4:6) = X0sampPostLTM(4:6) + deltaVLTM;
    [T_LTMtoLCM,X_LTMtoLCM] = ode45(@projectDynamics, LTM:step:LCM, ...
                                    X0sampPostLTM, options, params, sc);
    deltaXsampLCM = [X_LTMtoLCM(end,1:3) - Xtotal(indexLCM,1:3), ...
                     X_LTMtoLCM(end,4:6) - Xtotal(indexLCM,4:6)].';
    posMagError(i) = norm(deltaXsampLCM(1:3));
    velMagError(i) = norm(deltaXsampLCM(4:6));
    deltaVLCM = -inv(STM_LCM_Target(1:3,4:6))*STM_LCM_Target(1:3,1:3)* ...
                deltaXsampLCM(1:3) - deltaXsampLCM(4:6);
    deltaVLCMmag(i) = norm(deltaVLCM);
    X0sampPostLCM      = X_LTMtoLCM(end,:)';
    X0sampPostLCM(4:6) = X0sampPostLCM(4:6) + deltaVLCM;
    [TpostLCM,XpostLCM]  = ode45(@projectDynamics, ...
                                 LCM:step:LCM+227*params.tday, ...
                                 X0sampPostLCM, options, params, sc);
    Xsamp = [XpreLTM; X_LTMtoLCM(2:end,:); XpostLCM(2:end,:)];
    Tsamp = [TpreLTM; T_LTMtoLCM(2:end);   TpostLCM(2:end)];
    Rerror(i) = norm(Xsamp(indexLCM_th,1:3) - Xtotal(indexLCM_th,1:3));
    Verror(i) = norm(Xsamp(indexLCM_th,4:6) - Xtotal(indexLCM_th,4:6));
end
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
step    = 10;
X0preLTM = X0fit;
[TpreLTM, XpreLTM] = ode45(@projectDynamics, ET0:step:LTM, ...
                           X0preLTM, options, params, sc);
X0_LTM = XpreLTM(end,:).';
[Tpost_no, Xpost_no] = ode45(@projectDynamics, LTM:step:LTM+3600, ...
                             X0_LTM, options, params, sc);
XnoLTM = [XpreLTM(end-360:end,:); Xpost_no(2:end,:)];
TnoLTM = [TpreLTM(end-360:end);    Tpost_no(2:end)];
deltaVLTM = [0.6931075; -0.8462091; 0.0956979];
X0_LTM_burn        = X0_LTM;
X0_LTM_burn(4:6)   = X0_LTM_burn(4:6) + deltaVLTM;
[Tpost_yes, Xpost_yes] = ode45(@projectDynamics, LTM:step:LTM+3600, ...
                               X0_LTM_burn, options, params, sc);
XwithLTM = [XpreLTM(end-360:end,:); Xpost_yes(2:end,:)];
TwithLTM = [TpreLTM(end-360:end);    Tpost_yes(2:end)];
xEarth_LTM = cspice_spkezr('Earth', TnoLTM.', 'ECLIPJ2000','NONE','Sun');  % 6xN_noLTM
xEarth_LTM = xEarth_LTM.';  % N_noLTM x 6
XrelnoLTM   = XnoLTM(:,1:6)   - xEarth_LTM;
XrelwithLTM = XwithLTM(:,1:6) - xEarth_LTM;
R_EMOtoEME = inv(params.R_EMEtoEMO);
XrelnoLTMEME   = XrelnoLTM;
XrelnoLTMEME(:,1:3) = (R_EMOtoEME * XrelnoLTM(:,1:3).').';
XrelnoLTMEME(:,4:6) = (R_EMOtoEME * XrelnoLTM(:,4:6).').';
XrelwithLTMEME   = XrelwithLTM;
XrelwithLTMEME(:,1:3) = (R_EMOtoEME * XrelwithLTM(:,1:3).').';
XrelwithLTMEME(:,4:6) = (R_EMOtoEME * XrelwithLTM(:,4:6).').';
obsNoLTM = computeObservations(TnoLTM,   XrelnoLTMEME,   X0fit, params, st);
obsNoLTM(:,3) = obsNoLTM(:,3) + randn(size(obsNoLTM,1),1)*1e-3;
obsNoLTM(:,4) = obsNoLTM(:,4) + randn(size(obsNoLTM,1),1)*1e-7;
obsWithLTM = computeObservations(TwithLTM, XrelwithLTMEME, X0fit, params, st);
obsWithLTM(:,3) = obsWithLTM(:,3) + randn(size(obsWithLTM,1),1)*1e-3;
obsWithLTM(:,4) = obsWithLTM(:,4) + randn(size(obsWithLTM,1),1)*1e-7;
residuals = obsWithLTM - obsNoLTM;
optsV = struct();
r_sc_sun    = Xtotal(:,1:3);
r_earth_sun = xEarth(:,1:3);
r_moon_sun  = xEarth(:,1:3) + xMoon(:,1:3);
optsH = optsV;
optsH.videoFile   = 'mission_Heliocentric_followcam.mp4';
optsH.primaryName = 'Sun';
optsH.body2Name   = 'Earth';
optsH.showBody3   = true;
optsH.body3Name   = 'Moon';
optsH.followMode = 'axis';
optsH.r_body3     = r_moon_sun;
optsH.msBody         = 120;
optsH.showRangePanel = true;
writeMissionFollowCamVideo(Ttotal, r_sc_sun, r_earth_sun, ET0, LTM, LCM, deltaVLTM_nom, optsH);
function drawTexturedSphere(ax, center, R, texPath, res, fallbackColor)
[X,Y,Z] = sphere(max(10,res));
X = R*X + center(1); Y = R*Y + center(2); Z = R*Z + center(3);
h = surf(ax, X, Y, Z, 'EdgeColor','none', 'HandleVisibility','off');
if ~isempty(texPath) && exist(texPath,'file')
    img = flipud(imread(texPath));
    set(h, 'FaceColor','texturemap', 'CData', img);
else
    set(h, 'FaceColor', fallbackColor);
end
end
function opts = fillDefaultsFull(opts)
def.earthRadius_km = 6371;
def.moonRadius_km  = 1737;
def.earthTex = 'earth_texture.jpg';
def.moonTex  = 'moon_texture.jpg';
def.earthFallbackColor = [0.2 0.4 0.9];
def.moonFallbackColor  = [0.7 0.7 0.7];
def.sphereRes = 70;
def.maxPoints = 12000;
def.maxMoonPoints = 4000;
def.bodyFrac = 0.02;
def.arrowFrac = 0.03;
def.moonArcDays = 20;
def.trajLW = 1.7;
def.moonLW = 1.2;
def.viewVec = [-135 25];
def.lim = [];
def.limPad = 1.05;
opts = applyDefaults(opts, def);
end
function opts = fillDefaultsZoom(opts)
def.earthRadius_km = 6371;
def.moonRadius_km  = 1737;
def.earthTex = 'earth_texture.jpg';
def.moonTex  = 'moon_texture.jpg';
def.earthFallbackColor = [0.2 0.4 0.9];
def.moonFallbackColor  = [0.7 0.7 0.7];
def.sphereRes = 70;
def.maxPoints = 6000;
def.maxMoonPoints = 3000;
def.windowDays = 10;
def.moonArcDays = 27.3;
def.bodyFrac = 0.08;
def.arrowFrac = 0.05;
def.trajLW = 2.0;
def.moonLW = 1.4;
def.viewVec = [-135 25];
def.limPad = 1.10;
opts = applyDefaults(opts, def);
end
function opts = applyDefaults(opts, def)
fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(opts,f) || isempty(opts.(f))
        opts.(f) = def.(f);
    end
end
end