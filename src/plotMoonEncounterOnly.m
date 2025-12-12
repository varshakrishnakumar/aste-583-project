function [hFig, info] = plotMoonEncounterOnly(Ttotal, r_sc_E, r_moon_E, ET0, opts)
if nargin < 5 || isempty(opts), opts = struct(); end
def.windowDays       = 2;
def.focusRange_km    = [];
def.focusRangeFactor = 8;
def.focusRangeMin_km = 8*1737.4;
def.focusRangeMax_km = 2e5;
def.minSCPoints      = 80;
def.moonArcDays          = 1.5;
def.showMoonMotion       = true;
def.moonMotionMode       = 'scaled';
def.moonMotionRadiusFrac = 0.70;
def.moonLW               = 1.4;
def.moonTex          = 'moon_texture.jpg';
def.sphereRes        = 80;
def.moonScale        = 1.0;
def.moonRadius_km    = 1737.4;
def.maxMoonFracOfCA  = 0.85;
def.theme            = 'dark';
def.addStars         = true;
def.nStars           = 900;
def.colorByTime      = true;
def.showColorbar     = true;
def.cmap             = 'blueorange';
def.trajLW           = 3.0;
def.trajLWGlow       = 7.0;
def.moonLW           = 1.4;
def.viewVec          = [-120 20];
def.limPad           = 1.12;
def.showCALine       = true;
def.showVelArrow     = true;
def.arrowFrac        = 0.18;
def.fontName         = 'Helvetica';
def.fontSize         = 16;
def.gridAlpha        = 0.10;
def.titlePrefix      = 'Lunar Flyby Close-up';
def.useSpiceUtc      = true;
def.utcFormat        = 'ISOC';
def.utcPrec          = 3;
def.subtitleET0      = true;
def.camvaDeg = 8;
opts = applyDefaults(opts, def);
T = Ttotal(:);
N = numel(T);
r_sc_E   = r_sc_E(:,1:3);
r_moon_E = r_moon_E(:,1:3);
if size(r_sc_E,1) ~= N || size(r_moon_E,1) ~= N
    error('Ttotal, r_sc_E, and r_moon_E must have the same number of rows.');
end
if N < 3
    error('Need at least 3 time samples.');
end
tLooksLikeSeconds = max(abs(T)) > 1e5;
if tLooksLikeSeconds
    tSec = T;
else
    tSec = T * 86400;
end
daySec = 86400;
r_rel = r_sc_E - r_moon_E;
d = vecnorm(r_rel,2,2);
[dmin_grid, i0] = min(d);
[tCA_sec, dmin_km, caRefined] = refineCA_quadratic(tSec, d, i0);
r_rel_CA = interp1(tSec, r_rel, tCA_sec, 'pchip');
alt_km = dmin_km - opts.moonRadius_km;
useProvidedV = isfield(opts,'v_sc_E') && isfield(opts,'v_moon_E') && ...
               ~isempty(opts.v_sc_E) && ~isempty(opts.v_moon_E) && ...
               size(opts.v_sc_E,1)==N && size(opts.v_moon_E,1)==N;
if useProvidedV
    v_rel = opts.v_sc_E(:,1:3) - opts.v_moon_E(:,1:3);
else
    v_rel = centralDiff(tSec, r_rel);
end
v_rel_CA  = safeInterp(tSec, v_rel, tCA_sec);
v_rel_mag = norm(v_rel_CA);
v_rad = NaN; v_tan = NaN; phi_deg = NaN;
if isfinite(v_rel_mag) && v_rel_mag > 0
    rhat = r_rel_CA / max(norm(r_rel_CA), eps);
    v_rad = dot(v_rel_CA, rhat);
    v_tan = sqrt(max(v_rel_mag^2 - v_rad^2, 0));
    cosphi = dot(r_rel_CA, v_rel_CA) / (max(norm(r_rel_CA),eps)*max(v_rel_mag,eps));
    cosphi = max(min(cosphi,1),-1);
    phi_deg = acosd(cosphi);
end
halfWinSC = (opts.windowDays*daySec)/2;
if ~isempty(opts.focusRange_km)
    range_km = opts.focusRange_km;
else
    range_km = opts.focusRangeFactor * dmin_km;
end
range_km = max(range_km, opts.focusRangeMin_km);
range_km = min(range_km, opts.focusRangeMax_km);
[jL, jR] = contiguousAround(i0, tSec, tCA_sec, halfWinSC, d, range_km);
relax = 0;
while (jR-jL+1) < opts.minSCPoints && range_km < opts.focusRangeMax_km && relax < 8
    range_km = min(range_km*1.25, opts.focusRangeMax_km);
    [jL, jR] = contiguousAround(i0, tSec, tCA_sec, halfWinSC, d, range_km);
    relax = relax + 1;
end
idxSC = false(N,1); idxSC(jL:jR) = true;
if nnz(idxSC) < 2
    idxSC(max(i0-1,1):min(i0+1,N)) = true;
end
r_sc_seg = r_rel(idxSC,:);
t_seg    = tSec(idxSC);
dt_hr    = (t_seg - tCA_sec)/3600;
Rm_true = opts.moonRadius_km;
Rm_plot = Rm_true * opts.moonScale;
Rm_plot = min(Rm_plot, opts.maxMoonFracOfCA * dmin_km);
radSC = vecnorm(r_sc_seg,2,2);
radSC = radSC(isfinite(radSC));
if isempty(radSC), radSC = dmin_km; end
limData = robustPercentile(radSC, 99.5);
lim = max([1.15*Rm_plot, limData]);
lim = lim * opts.limPad;
r_moon_hint = [];
if opts.showMoonMotion
    r_moon_CA_E = interp1(tSec, r_moon_E, tCA_sec, 'pchip');
    halfWinMoon = (opts.moonArcDays*daySec)/2;
    idxHint = abs(tSec - tCA_sec) <= halfWinMoon;
    dispE = r_moon_E(idxHint,:) - r_moon_CA_E;
    if size(dispE,1) >= 2
        switch lower(string(opts.moonMotionMode))
            case "true"
                r_moon_hint = dispE;
            otherwise
                R_hint = opts.moonMotionRadiusFrac * lim;
                mags = vecnorm(dispE,2,2);
                mmax = max(mags);
                if isfinite(mmax) && mmax > 0
                    r_moon_hint = dispE * (R_hint / mmax);
                end
        end
    end
end
utcStr = '';
if opts.useSpiceUtc && tLooksLikeSeconds && exist('cspice_et2utc','file')==2
    try
        utcStr = cspice_et2utc(tCA_sec, opts.utcFormat, opts.utcPrec);
    catch
        utcStr = '';
    end
end
dtFromET0_days = [];
if ~isempty(ET0) && isscalar(ET0)
    if tLooksLikeSeconds
        dtFromET0_days = (tCA_sec - ET0)/86400;
    else
        dtFromET0_days = (tCA_sec/86400 - ET0);
    end
end
isDark = strcmpi(opts.theme,'dark');
if isDark
    bg = [0.03 0.03 0.05];
    fg = [0.95 0.95 0.97];
    gridC = [1 1 1]*0.18;
    trajGlowC = [0.9 0.9 1.0];
    moonArcC = [1 1 1]*0.35;
    caC = [1.0 0.45 0.15];
else
    bg = [1 1 1];
    fg = [0 0 0];
    gridC = [0 0 0]*0.25;
    trajGlowC = [0.2 0.2 0.2];
    moonArcC = [0.35 0.35 0.35];
    caC = [0.85 0.325 0.098];
end
hFig = figure('Color',bg,'Renderer','opengl');
set(hFig,'InvertHardcopy','off');
ax = axes('Parent',hFig); hold(ax,'on');
axis(ax,'equal'); axis(ax,'vis3d');
ax.Color = bg;
ax.FontName = opts.fontName;
ax.FontSize = opts.fontSize;
ax.XColor = fg; ax.YColor = fg; ax.ZColor = fg;
ax.GridColor = gridC; ax.MinorGridColor = gridC;
ax.GridAlpha = opts.gridAlpha; ax.MinorGridAlpha = opts.gridAlpha;
grid(ax,'on');
xlabel(ax,'X (km)  [Moon-centered]');
ylabel(ax,'Y (km)  [Moon-centered]');
zlabel(ax,'Z (km)  [Moon-centered]');
if isDark && opts.addStars
    addStarField(ax, lim, opts.nStars, fg);
end
line1 = opts.titlePrefix;
if ~isempty(utcStr)
    timePart = sprintf('CA: %s', utcStr);
else
    timePart = sprintf('CA @ t = %.6f days', tCA_sec/86400);
end
tPlus = '';
if ~isempty(dtFromET0_days) && opts.subtitleET0
    tPlus = sprintf(' (T+%.3f days)', dtFromET0_days);
end
line2 = sprintf('%s%s  |  d_{min}=%.1f km  |  alt=%.1f km', timePart, tPlus, dmin_km, alt_km);
if isfinite(v_rel_mag)
    line2 = sprintf('%s  |  v_{rel}=%.3f km/s', line2, v_rel_mag);
end
title(ax,{line1,line2},'FontWeight','bold','Color',fg);
drawTexturedSphere(ax, [0 0 0], Rm_plot, opts.moonTex, opts.sphereRes, isDark);
drawHalo(ax, [0 0 0], 1.02*Rm_plot, isDark);
hMoonHint = [];
if ~isempty(r_moon_hint)
    hMoonHint = plot3(ax, r_moon_hint(:,1), r_moon_hint(:,2), r_moon_hint(:,3), '--', ...
        'LineWidth', opts.moonLW, 'Color', moonArcC, 'DisplayName','Moon motion (hint)');
end
hTraj = [];
if opts.colorByTime && size(r_sc_seg,1) >= 2
    hGlow = plot3(ax, r_sc_seg(:,1), r_sc_seg(:,2), r_sc_seg(:,3), '-', ...
        'LineWidth', opts.trajLWGlow, 'Color', trajGlowC, 'HandleVisibility','off');
    x = r_sc_seg(:,1); y = r_sc_seg(:,2); z = r_sc_seg(:,3); c = dt_hr(:);
    hTraj = surface(ax, [x x], [y y], [z z], [c c], ...
        'FaceColor','none','EdgeColor','interp','LineWidth', opts.trajLW);
    hTraj.DisplayName = 'Spacecraft (time-colored)';
    dtMax = max(abs(dt_hr));
    if ~isfinite(dtMax) || dtMax < 1e-6, dtMax = 1; end
    caxis(ax, [-dtMax dtMax]);
    colormap(ax, chooseCmap(opts.cmap, 256, isDark));
    if opts.showColorbar
        cb = colorbar(ax);
        cb.Color = fg;
        cb.Label.String = 'Hours from CA';
        cb.Label.Color  = fg;
    end
else
    hTraj = plot3(ax, r_sc_seg(:,1), r_sc_seg(:,2), r_sc_seg(:,3), '-', ...
        'LineWidth', opts.trajLW, 'Color', trajGlowC, 'DisplayName','Spacecraft');
end
posCA = r_rel_CA;
hCA = scatter3(ax, posCA(1), posCA(2), posCA(3), 80, caC, 'filled', ...
    'MarkerEdgeColor', fg, 'LineWidth', 1.0, 'DisplayName','Closest approach');
caLineC = 0.55*fg + 0.45*bg;
plot3(ax, [0 posCA(1)], [0 posCA(2)], [0 posCA(3)], ':', ...
    'Color', caLineC, 'LineWidth', 1.2, 'HandleVisibility','off');
if opts.showVelArrow && isfinite(v_rel_mag) && v_rel_mag > 0
    u = v_rel_CA / v_rel_mag;
    L = opts.arrowFrac * lim;
    quiver3(ax, posCA(1),posCA(2),posCA(3), L*u(1),L*u(2),L*u(3), 0, ...
        'Color', caC, 'LineWidth', 2.0, 'MaxHeadSize', 0.6, ...
        'HandleVisibility','off');
end
off = 0.06*lim;
txtBg = isDark*[0 0 0] + (~isDark)*[1 1 1];
text(ax, posCA(1)+off, posCA(2)+off, posCA(3)+off, ...
    sprintf('CA: %.1f km', dmin_km), ...
    'Color', fg, 'FontWeight','bold', ...
    'BackgroundColor', txtBg, 'Margin', 4, 'Clipping','on');
xlim(ax, [-lim lim]); ylim(ax, [-lim lim]); zlim(ax, [-lim lim]);
view(ax, opts.viewVec);
camtarget(ax, [0 0 0]);
camva(ax, opts.camvaDeg);
lighting(ax,'gouraud'); material(ax,'dull');
camlight(ax,'headlight');
camlight(ax,'right');
legItems = [hTraj, hCA];
if ~isempty(hMoonHint), legItems = [hMoonHint, legItems]; end
lg = legend(ax, legItems, 'Location','northeast');
lg.TextColor = fg;
lg.Color = txtBg;
lg.EdgeColor = gridC;
lines = {};
if ~isempty(utcStr), lines{end+1} = sprintf('UTC CA: %s', utcStr); end
if ~isempty(dtFromET0_days), lines{end+1} = sprintf('T+%.6f days', dtFromET0_days); end
lines{end+1} = sprintf('dmin: %.3f km', dmin_km);
lines{end+1} = sprintf('alt:  %.3f km', alt_km);
if isfinite(v_rel_mag)
    lines{end+1} = sprintf('v_rel: %.4f km/s', v_rel_mag);
    lines{end+1} = sprintf('v_rad: %.4f km/s', v_rad);
    lines{end+1} = sprintf('v_tan: %.4f km/s', v_tan);
    lines{end+1} = sprintf('angle(r,v): %.2f deg', phi_deg);
end
lines{end+1} = ternary(caRefined, 'CA time refined (quadratic)', 'CA time on grid');
annotation(hFig,'textbox', [0.02 0.64 0.34 0.28], ...
    'String', lines, 'FitBoxToText','on', ...
    'BackgroundColor', txtBg, 'EdgeColor', gridC, ...
    'Color', fg, 'FontName', opts.fontName, ...
    'FontSize', max(opts.fontSize-1, 9), 'Interpreter','none');
info = struct();
info.tCA_sec      = tCA_sec;
info.tCA_days     = tCA_sec/86400;
info.utc          = utcStr;
info.dmin_km      = dmin_km;
info.dmin_grid_km = dmin_grid;
info.alt_km       = alt_km;
info.v_rel_CA     = v_rel_CA;
info.v_rel_kms    = v_rel_mag;
info.v_rad_kms    = v_rad;
info.v_tan_kms    = v_tan;
info.angle_rv_deg = phi_deg;
info.caRefined    = caRefined;
info.idxCA_grid   = i0;
info.focusRange_km = range_km;
info.segmentIdx   = [jL jR];
end
function [jL, jR] = contiguousAround(i0, tSec, tCA, halfWin, dist, maxDist)
jL = i0; jR = i0;
while jL > 1
    if abs(tSec(jL-1)-tCA) <= halfWin && dist(jL-1) <= maxDist
        jL = jL - 1;
    else
        break
    end
end
while jR < numel(tSec)
    if abs(tSec(jR+1)-tCA) <= halfWin && dist(jR+1) <= maxDist
        jR = jR + 1;
    else
        break
    end
end
end
function vq = safeInterp(t, v, tq)
try
    vq = interp1(t, v, tq, 'pchip');
catch
    vq = [NaN NaN NaN];
end
end
function drawTexturedSphere(axh, center, R, texPath, res, isDark)
[X,Y,Z] = sphere(max(16,res));
X = R*X + center(1); Y = R*Y + center(2); Z = R*Z + center(3);
h = surf(axh, X, Y, Z, 'EdgeColor','none', 'HandleVisibility','off');
h.AmbientStrength  = 0.45;
h.DiffuseStrength  = 0.70;
h.SpecularStrength = 0.10;
if ~isempty(texPath) && exist(texPath,'file')
    img = flipud(imread(texPath));
    set(h, 'FaceColor','texturemap', 'CData', img);
else
    set(h, 'FaceColor', ternary(isDark,[0.7 0.7 0.75],[0.75 0.75 0.75]));
end
end
function drawHalo(axh, center, R, isDark)
[X,Y,Z] = sphere(40);
X = R*X + center(1); Y = R*Y + center(2); Z = R*Z + center(3);
hc = surf(axh, X, Y, Z, 'EdgeColor','none', 'HandleVisibility','off');
hc.FaceColor = ternary(isDark,[0.7 0.8 1.0],[0.9 0.9 1.0]);
hc.FaceAlpha = ternary(isDark,0.08,0.04);
end
function addStarField(ax, lim, nStars, fg)
M = max(3*nStars, 2000);
P = (rand(M,3)*2 - 1) * lim;
r = vecnorm(P,2,2);
P = P(r > 0.78*lim,:);
if size(P,1) > nStars
    P = P(1:nStars,:);
end
s = 4 + 6*rand(size(P,1),1);
h = scatter3(ax, P(:,1),P(:,2),P(:,3), s, fg, 'filled', 'HandleVisibility','off');
try
    h.MarkerFaceAlpha = 0.18;
    h.MarkerEdgeAlpha = 0.0;
catch
end
end
function cmap = chooseCmap(name, n, isDark)
if nargin < 2, n = 256; end
if nargin < 3, isDark = false; end
name = lower(string(name));
switch name
    case "blueorange"
        cmap = blueOrange(n, isDark);
    case "turbo"
        if exist('turbo','file')==2 || exist('turbo','builtin')==5
            cmap = turbo(n);
        else
            cmap = parula(n);
        end
    case "parula"
        cmap = parula(n);
    otherwise
        cmap = parula(n);
end
end
function cmap = blueOrange(n, isDark)
blue   = [0.20 0.55 0.95];
orange = [1.00 0.55 0.18];
mid    = ternary(isDark, [0.92 0.92 0.96], [0.15 0.15 0.18]);
n1 = floor(n/2);
n2 = n - n1;
a = lerpRow(blue, mid, n1);
b = lerpRow(mid, orange, n2);
cmap = [a; b];
end
function A = lerpRow(c1, c2, n)
t = linspace(0,1,max(n,2)).';
A = c1.*(1-t) + c2.*t;
A = A(1:n,:);
end
function [tCA_sec, dmin_km, refined] = refineCA_quadratic(tSec, d, i0)
N = numel(tSec);
iFit = max(1,i0-2) : min(N,i0+2);
refined = false;
t0 = tSec(i0);
if numel(iFit) < 3
    tCA_sec = t0;
    dmin_km = d(i0);
    return
end
dt = tSec(iFit) - t0;
s  = d(iFit).^2;
p = polyfit(dt, s, 2);
a = p(1); b = p(2);
if ~(isfinite(a) && a > 0)
    tCA_sec = t0; dmin_km = d(i0);
    return
end
dtMin = -b/(2*a);
if dtMin < min(dt) || dtMin > max(dt)
    tCA_sec = t0; dmin_km = d(i0);
    return
end
sMin = polyval(p, dtMin);
tCA_sec = t0 + dtMin;
dmin_km = sqrt(max(sMin,0));
refined = true;
end
function dXdt = centralDiff(t, X)
N = numel(t);
dXdt = zeros(size(X));
dXdt(1,:)   = (X(2,:)   - X(1,:))     / (t(2)   - t(1));
dXdt(end,:) = (X(end,:) - X(end-1,:)) / (t(end) - t(end-1));
for i = 2:N-1
    dXdt(i,:) = (X(i+1,:) - X(i-1,:)) / (t(i+1) - t(i-1));
end
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
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
function p = robustPercentile(x, pct)
x = x(:);
x = x(isfinite(x));
if isempty(x), p = NaN; return; end
x = sort(x);
k = 1 + (numel(x)-1) * (pct/100);
f = floor(k); c = ceil(k);
if f == c
    p = x(f);
else
    p = x(f) + (k-f) * (x(c)-x(f));
end
end