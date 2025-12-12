function hFig = plotLunarFlybyZoom(Ttotal, r_sc, r_moon, ET0, LTM, LCM, deltaVLTM, opts)
if nargin < 8 || isempty(opts), opts = struct(); end
opts = fillDefaultsZoom(opts);
T = Ttotal(:);
r_sc   = r_sc(:,1:3);
r_moon = r_moon(:,1:3);
dSM = vecnorm(r_sc - r_moon, 2, 2);
[dmin, iCA] = min(dSM);
tCA = T(iCA);
win = opts.windowDays*86400;
maskWin = abs(T - tCA) <= win/2;
idxWin = find(maskWin);
arc = opts.moonArcDays*86400;
maskMoon = abs(T - tCA) <= arc/2;
idxMoon = find(maskMoon);
rLocal = [r_sc(idxWin,:); r_moon(idxMoon,:); 0 0 0];
lim = max(vecnorm(rLocal,2,2));
lim = lim * opts.limPad;
scale = max(1, (opts.bodyFrac*lim)/opts.earthRadius_km);
Re = opts.earthRadius_km * scale;
Rm = opts.moonRadius_km  * scale;
stride = max(1, ceil(numel(idxWin)/opts.maxPoints));
idxWin = idxWin(1:stride:end);
strideM = max(1, ceil(numel(idxMoon)/opts.maxMoonPoints));
idxMoon = idxMoon(1:strideM:end);
hFig = figure('Color','w');
ax = axes('Parent',hFig); hold(ax,'on'); grid(ax,'on');
axis(ax,'equal'); axis(ax,'vis3d');
xlabel(ax,'X_E (km)'); ylabel(ax,'Y_E (km)'); zlabel(ax,'Z_E (km)');
title(ax, sprintf('Zoomed Flyby (centered on closest approach) | d_{min}=%.0f km @ %.2f days', ...
    dmin, (tCA-ET0)/86400));
drawTexturedSphere(ax, [0 0 0], Re, opts.earthTex, opts.sphereRes, opts.earthFallbackColor);
drawTexturedSphere(ax, r_moon(iCA,:), Rm, opts.moonTex, opts.sphereRes-15, opts.moonFallbackColor);
hMoonArc = plot3(ax, r_moon(idxMoon,1), r_moon(idxMoon,2), r_moon(idxMoon,3), ...
    '--', 'LineWidth', opts.moonLW, 'DisplayName','Moon (arc)');
hTraj = plot3(ax, r_sc(idxWin,1), r_sc(idxWin,2), r_sc(idxWin,3), ...
    '-', 'LineWidth', opts.trajLW, 'DisplayName','Spacecraft (window)');
hCA = plot3(ax, r_sc(iCA,1), r_sc(iCA,2), r_sc(iCA,3), ...
    'ko','MarkerSize',8,'LineWidth',2,'DisplayName','Closest approach');
posCA = r_sc(iCA,:);
off   = 0.03*lim;
text(ax, posCA(1)+off, posCA(2)+off, posCA(3)+off, ...
     sprintf('CA = %.0f km', dmin), ...
     'Color','k','FontWeight','bold','BackgroundColor','w','Margin',2);
plot3(ax, [r_sc(iCA,1) r_moon(iCA,1)], ...
         [r_sc(iCA,2) r_moon(iCA,2)], ...
         [r_sc(iCA,3) r_moon(iCA,3)], ...
      'k:', 'LineWidth', 1.5, 'HandleVisibility','off');
fprintf('Closest approach = %.1f km, Moon radius = %.1f km\n', dmin, opts.moonRadius_km);
[~, iLTM] = min(abs(T - LTM));
[~, iLCM] = min(abs(T - LCM));
if abs(T(iLTM) - tCA) <= win/2
    plot3(ax, r_sc(iLTM,1), r_sc(iLTM,2), r_sc(iLTM,3), 'ro', 'MarkerSize',8,'LineWidth',2,'DisplayName','LTM');
    dv = deltaVLTM(:).';
    dv_hat = dv / max(norm(dv), eps);
    arrowLen = opts.arrowFrac * lim;
    quiver3(ax, r_sc(iLTM,1), r_sc(iLTM,2), r_sc(iLTM,3), ...
        dv_hat(1)*arrowLen, dv_hat(2)*arrowLen, dv_hat(3)*arrowLen, 0, ...
        'LineWidth',2,'MaxHeadSize',0.6,'DisplayName','\DeltaV @ LTM');
end
if abs(T(iLCM) - tCA) <= win/2
    plot3(ax, r_sc(iLCM,1), r_sc(iLCM,2), r_sc(iLCM,3), 'go', 'MarkerSize',8,'LineWidth',2,'DisplayName','LCM');
end
mins = min([r_sc(idxWin,:); r_moon(idxMoon,:); 0 0 0], [], 1);
maxs = max([r_sc(idxWin,:); r_moon(idxMoon,:); 0 0 0], [], 1);
pad  = 0.08*(maxs - mins);
xlim(ax, [mins(1)-pad(1), maxs(1)+pad(1)]);
ylim(ax, [mins(2)-pad(2), maxs(2)+pad(2)]);
zlim(ax, [mins(3)-pad(3), maxs(3)+pad(3)]);
view(ax, opts.viewVec);
lighting(ax,'gouraud'); material(ax,'dull'); camlight(ax,'headlight');
legend(ax,'Location','northeast');
end
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
def.moonArcDays          = 1.5;
def.showMoonMotion       = true;
def.moonMotionMode       = 'scaled';
def.moonMotionRadiusFrac = 0.70;
def.moonLW               = 1.4;
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
def.bodyFrac = 0.08;
def.arrowFrac = 0.05;
def.moonArcDays          = 1.5;
def.showMoonMotion       = true;
def.moonMotionMode       = 'scaled';
def.moonMotionRadiusFrac = 0.70;
def.moonLW               = 1.4;
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