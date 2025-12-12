function [hFig, info] = plotMission3DConops(Ttotal, r_sc_E, r_moon_E, ET0, LTM, LCM, deltaVLTM, opts)
if nargin < 8 || isempty(opts), opts = struct(); end
Ttotal = Ttotal(:);
r_sc_E = r_sc_E(:,1:3);
def.frame          = 'heliocentric';
def.units          = 'AU';
def.rEarthSun      = [];
def.useSpiceEarth  = false;
def.spiceFrame     = 'ECLIPJ2000';
def.spiceCenter    = 'Sun';
def.spiceCorr      = 'NONE';
def.scIsHelio      = false;
def.theme          = 'dark';
def.addStars       = true;
def.showRangePanel = true;
def.showEarthTrack = true;
def.showSun        = true;
def.show1AUring    = true;
def.maxPoints      = 12000;
def.maxEarthPoints = 4000;
def.centerMode     = 'sun';
def.limPad         = 1.08;
def.viewVec        = [-135 25];
def.fontName       = 'Helvetica';
def.fontSize       = 11;
def.lwFaint   = 1.2;
def.lwMain    = 2.8;
def.cPre      = [0.65 0.65 0.65];
def.cMid      = [0.95 0.45 0.10];
def.cPost     = [0.20 0.70 1.00];
def.msEvent   = 8;
def.msBody    = 70;
def.cSun      = [1.00 0.80 0.10];
def.cEarth    = [0.20 0.60 1.00];
def.earthTrackStyle = ':';
def.earthTrackLW    = 1.1;
def.dvFrac       = 0.035;
def.dvColor      = [0.75 0.10 0.95];
def.dvLW         = 2.4;
def.dvHead       = 1.2;
def.dvOffsetFrac = 0.010;
def.useSpiceUtc   = true;
def.utcFormat     = 'ISOC';
def.utcPrec       = 0;
def.titlePrefix   = 'Heliocentric Mission CONOPS';
def.useOrtho      = true;
def.showInset       = true;
def.insetDaysPad    = 6;
def.insetSizeFrac   = 0.36;
def.insetCorner     = 'northeast';
def.insetPad        = 1.25;
def.maxInsetPoints  = 3500;
def.insetShowEarth  = true;
def.insetTitle      = 'Zoom: LTM / LCM';
def.showMetricsBox  = true;
opts = applyDefaults(opts, def);
T = Ttotal;
tLooksLikeSeconds = max(abs(T)) > 1e5;
if tLooksLikeSeconds
    tSec = T;
    ET0_sec = ET0;
    LTM_sec = LTM;
    LCM_sec = LCM;
    tDays = (tSec - ET0_sec)/86400;
else
    tSec = T * 86400;
    ET0_sec = ET0 * 86400;
    LTM_sec = LTM * 86400;
    LCM_sec = LCM * 86400;
    tDays = (T - ET0);
end
N = size(r_sc_E,1);
iLTM = nearestIndex(tSec, LTM_sec);
iLCM = nearestIndex(tSec, LCM_sec);
r_earth_sun = zeros(N,3);
earthSrc = "none";
if strcmpi(opts.frame,'heliocentric') && ~opts.scIsHelio
    if ~isempty(opts.rEarthSun)
        r_earth_sun = opts.rEarthSun(:,1:3);
        earthSrc = "opts.rEarthSun";
        if size(r_earth_sun,1) ~= N
            error('opts.rEarthSun must be Nx3 to match Ttotal.');
        end
    elseif opts.useSpiceEarth
        if ~(tLooksLikeSeconds && exist('cspice_spkezr','file')==2)
            error('useSpiceEarth requires ET seconds in Ttotal and cspice_spkezr on path.');
        end
        xE = cspice_spkezr('Earth', tSec.', opts.spiceFrame, opts.spiceCorr, opts.spiceCenter); % 6xN
        r_earth_sun = xE(1:3,:).';
        earthSrc = "SPICE";
    else
        error(['Heliocentric frame needs Earth wrt Sun. Provide opts.rEarthSun (Nx3 km) ' ...
               'or set opts.useSpiceEarth=true (ET seconds + SPICE).']);
    end
end
AU_km = 149597870.7;
switch lower(char(string(opts.frame)))
    case 'heliocentric'
        if opts.scIsHelio
            r_sc_sun = r_sc_E;
        else
            r_sc_sun = r_earth_sun + r_sc_E;
        end
        r_plot_sc    = r_sc_sun;
        r_plot_earth = r_earth_sun;
        originPos    = [0 0 0];
        originName   = 'Sun';
        if strcmpi(opts.units,'au')
            scale = 1/AU_km;
            unitLabel = 'AU';
        else
            scale = 1;
            unitLabel = 'km';
        end
    otherwise
        r_sc_sun      = nan(N,3);
        r_plot_sc     = r_sc_E;
        r_plot_earth  = zeros(N,3);
        originPos     = [0 0 0];
        originName    = 'Earth';
        scale = 1;
        unitLabel = 'km';
end
r_plot_sc    = r_plot_sc * scale;
r_plot_earth = r_plot_earth * scale;
originPos    = originPos * scale;
dv = deltaVLTM(:);
dvMag_mps = 1000*norm(dv);
rE_km = vecnorm(r_sc_E,2,2);
if strcmpi(opts.frame,'heliocentric')
    rS_AU = vecnorm(r_sc_sun,2,2)/AU_km;
else
    rS_AU = nan(N,1);
end
tLTMd = tDays(iLTM);
tLCMd = tDays(iLCM);
utcLTM = ''; utcLCM = '';
if opts.useSpiceUtc && tLooksLikeSeconds && exist('cspice_et2utc','file')==2
    try, utcLTM = cspice_et2utc(LTM_sec, opts.utcFormat, opts.utcPrec); end
    try, utcLCM = cspice_et2utc(LCM_sec, opts.utcFormat, opts.utcPrec); end
end
idxSC    = decimateWithAnchors(N, opts.maxPoints,      [1 iLTM iLCM N]);
idxEarth = decimateWithAnchors(N, opts.maxEarthPoints, [1 iLTM iLCM N]);
idxPre  = unique([idxSC(idxSC <= iLTM); iLTM]);
idxMid  = unique([idxSC(idxSC >= iLTM & idxSC <= iLCM); iLTM; iLCM]);
idxPost = unique([idxSC(idxSC >= iLCM); iLCM]);
[hFig, ax3, axR] = makeLayout(opts);
[bg, fg, gridC] = themeColors(opts.theme);
set(hFig,'Color',bg,'InvertHardcopy','off');
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
axis(ax3,'equal'); axis(ax3,'vis3d');
view(ax3, opts.viewVec);
set(ax3,'Color',bg,'XColor',fg,'YColor',fg,'ZColor',fg, ...
    'GridColor',gridC,'GridAlpha',0.25, ...
    'FontName',opts.fontName,'FontSize',opts.fontSize);
xlabel(ax3, sprintf('X (%s)', unitLabel));
ylabel(ax3, sprintf('Y (%s)', unitLabel));
zlabel(ax3, sprintf('Z (%s)', unitLabel));
if opts.useOrtho
    try, camproj(ax3,'orthographic'); end
    try, ax3.Projection = 'orthographic'; end
end
if isprop(ax3,'ClippingStyle')
    ax3.ClippingStyle = 'rectangle';
end
if strcmpi(opts.theme,'dark') && opts.addStars
    addStarField(ax3, r_plot_sc(idxSC,:), 800, fg);
end
h1AU = [];
if strcmpi(opts.frame,'heliocentric') && strcmpi(opts.units,'au') && opts.show1AUring
    th = linspace(0,2*pi,400);
    h1AU = plot3(ax3, cos(th), sin(th), 0*th, ':', 'Color', gridC, 'LineWidth', 1.1, ...
        'DisplayName','1 AU ring');
end
hOrigin = [];
if strcmpi(opts.frame,'heliocentric') && opts.showSun
    hOrigin = scatter3(ax3, originPos(1),originPos(2),originPos(3), opts.msBody, 'filled', ...
        'MarkerFaceColor', opts.cSun, 'MarkerEdgeColor', fg, 'DisplayName',originName);
    text(ax3, 0,0,0, ['  ' originName], 'Color', fg, 'FontWeight','bold');
end
hEarthTrack = [];
if strcmpi(opts.frame,'heliocentric') && opts.showEarthTrack
    hEarthTrack = plot3(ax3, r_plot_earth(idxEarth,1), r_plot_earth(idxEarth,2), r_plot_earth(idxEarth,3), ...
        opts.earthTrackStyle, 'Color', gridC, 'LineWidth', opts.earthTrackLW, 'DisplayName','Earth orbit');
    scatter3(ax3, r_plot_earth(iLTM,1), r_plot_earth(iLTM,2), r_plot_earth(iLTM,3), ...
        45, opts.cEarth, 'filled', 'MarkerEdgeColor', fg, 'HandleVisibility','off');
    text(ax3, r_plot_earth(iLTM,1), r_plot_earth(iLTM,2), r_plot_earth(iLTM,3), '  Earth @ LTM', ...
        'Color', fg, 'FontWeight','bold');
end
plot3(ax3, r_plot_sc(idxSC,1), r_plot_sc(idxSC,2), r_plot_sc(idxSC,3), '-', ...
    'Color', gridC, 'LineWidth', 1.0, 'HandleVisibility','off');
hPre = plot3(ax3, r_plot_sc(idxPre,1), r_plot_sc(idxPre,2), r_plot_sc(idxPre,3), '-', ...
    'Color', opts.cPre, 'LineWidth', opts.lwFaint, 'DisplayName','ET0 \rightarrow LTM');
hMid = plot3(ax3, r_plot_sc(idxMid,1), r_plot_sc(idxMid,2), r_plot_sc(idxMid,3), '-', ...
    'Color', opts.cMid, 'LineWidth', opts.lwMain, 'DisplayName','LTM \rightarrow LCM');
hPost = plot3(ax3, r_plot_sc(idxPost,1), r_plot_sc(idxPost,2), r_plot_sc(idxPost,3), '-', ...
    'Color', opts.cPost, 'LineWidth', opts.lwMain, 'DisplayName','LCM \rightarrow End');
hLTM = plot3(ax3, r_plot_sc(iLTM,1), r_plot_sc(iLTM,2), r_plot_sc(iLTM,3), 'o', ...
    'MarkerSize', opts.msEvent, 'MarkerFaceColor', [1 0.2 0.2], ...
    'MarkerEdgeColor', fg, 'LineWidth', 1.2, 'DisplayName','LTM');
hLCM = plot3(ax3, r_plot_sc(iLCM,1), r_plot_sc(iLCM,2), r_plot_sc(iLCM,3), 'o', ...
    'MarkerSize', opts.msEvent, 'MarkerFaceColor', [0.2 1.0 0.2], ...
    'MarkerEdgeColor', fg, 'LineWidth', 1.2, 'DisplayName','LCM');
text(ax3, r_plot_sc(iLTM,1), r_plot_sc(iLTM,2), r_plot_sc(iLTM,3), '  LTM', 'Color', fg, 'FontWeight','bold');
text(ax3, r_plot_sc(iLCM,1), r_plot_sc(iLCM,2), r_plot_sc(iLCM,3), '  LCM', 'Color', fg, 'FontWeight','bold');
ptsForBox = [r_plot_sc(idxSC,:); originPos];
if ~isempty(hEarthTrack), ptsForBox = [ptsForBox; r_plot_earth(idxEarth,:)]; end
[xl, yl, zl, boxSize] = computeBoxLimits(ptsForBox, opts);
xlim(ax3, xl); ylim(ax3, yl); zlim(ax3, zl);
dvHat = dv(:).';
if any(~isfinite(dvHat)) || norm(dvHat) < eps
    dvHat = [1 0 0];
else
    dvHat = dvHat / norm(dvHat);
end
arrowLen = opts.dvFrac * boxSize;
r0 = r_plot_sc(iLTM,:);
if norm(r0) > eps, rhat = r0 / norm(r0); else, rhat = dvHat; end
base = r0 + rhat * (opts.dvOffsetFrac * boxSize);
hDV = quiver3(ax3, base(1), base(2), base(3), ...
    dvHat(1)*arrowLen, dvHat(2)*arrowLen, dvHat(3)*arrowLen, 0, ...
    'Color', opts.dvColor, 'LineWidth', opts.dvLW, ...
    'MaxHeadSize', opts.dvHead, 'DisplayName','\DeltaV @ LTM');
line1 = opts.titlePrefix;
if ~isempty(utcLTM) && ~isempty(utcLCM)
    line2 = sprintf('LTM: %s  |  LCM: %s  |  \\DeltaV=%.0f m/s', utcLTM, utcLCM, dvMag_mps);
else
    line2 = sprintf('LTM: T+%.2f d  |  LCM: T+%.2f d  |  \\DeltaV=%.0f m/s', tLTMd, tLCMd, dvMag_mps);
end
title(ax3, {line1, line2}, 'Color', fg, 'FontWeight','bold');
legItems = [hPre hMid hPost hLTM hLCM hDV];
if ~isempty(hEarthTrack), legItems = [hEarthTrack legItems]; end
if ~isempty(h1AU), legItems = [h1AU legItems]; end
leg = legend(ax3, legItems, 'Location','northeastoutside');
leg.TextColor = fg; leg.EdgeColor = gridC; leg.Color = bg;
if opts.showInset
    axInset = addInsetLTM_LCM(hFig, ax3, opts, bg, fg, gridC, ...
        tSec, LTM_sec, LCM_sec, ...
        r_plot_sc, r_plot_earth, ...
        iLTM, iLCM, dvHat, unitLabel);
    if opts.useOrtho
        try, camproj(axInset,'orthographic'); end
        try, axInset.Projection = 'orthographic'; end
    end
    if isprop(axInset,'ClippingStyle')
        axInset.ClippingStyle = 'rectangle';
    end
end
if opts.showMetricsBox && strcmpi(opts.frame,'heliocentric')
    s = sprintf([ ...
        'Key metrics\n' ...
        'Earth source: %s | Units: %s\n' ...
        'LTM:  t=%.2f d | |r_{SC-E}|=%.2f Mm | |r_{SC-S}|=%.4f AU\n' ...
        'LCM:  t=%.2f d | |r_{SC-E}|=%.2f Mm | |r_{SC-S}|=%.4f AU\n' ...
        '\\DeltaV@LTM: %.0f m/s | TOF(LTM\\rightarrowLCM)=%.2f d'], ...
        char(earthSrc), unitLabel, ...
        tLTMd, rE_km(iLTM)/1e6, rS_AU(iLTM), ...
        tLCMd, rE_km(iLCM)/1e6, rS_AU(iLCM), ...
        dvMag_mps, (tLCMd - tLTMd));
    annotation(hFig,'textbox', [0.02 0.70 0.46 0.23], 'String', s, ...
        'FitBoxToText','on', 'BackgroundColor', bg, 'EdgeColor', gridC, ...
        'Color', fg, 'FontName', opts.fontName, 'FontSize', opts.fontSize-1);
end
if opts.showRangePanel && ~isempty(axR) && isvalid(axR)
    hold(axR,'on'); grid(axR,'on'); box(axR,'on');
    set(axR,'Color',bg,'XColor',fg,'YColor',fg, ...
        'GridColor',gridC,'GridAlpha',0.25, ...
        'FontName',opts.fontName,'FontSize',opts.fontSize);
    yyaxis(axR,'left');
    hRad = plot(axR, tDays, rS_AU, '-', 'LineWidth', 1.6, 'Color', [0.20 0.70 1.00]);
    ylabel(axR,'Heliocentric radius (AU)');
    axR.YAxis(1).Color = [0.20 0.70 1.00];
    yyaxis(axR,'right');
    hRE  = plot(axR, tDays, rE_km/1e6, '--', 'LineWidth', 1.6, 'Color', [1.00 0.55 0.10]);
    ylabel(axR,'Earth range (million km)');
    axR.YAxis(2).Color = [1.00 0.55 0.10];
    xline(axR, tLTMd, ':', 'LTM', 'Color', [1 0.2 0.2], 'LineWidth', 1.4);
    xline(axR, tLCMd, ':', 'LCM', 'Color', [0.2 0.9 0.2], 'LineWidth', 1.4);
    xlabel(axR,'Days since ET0');
    set([hRad hRE], 'Clipping','on');
end
info = struct();
info.frame = opts.frame;
info.units = unitLabel;
info.iLTM  = iLTM;
info.iLCM  = iLCM;
info.dv_mps = dvMag_mps;
end
function axI = addInsetLTM_LCM(hFig, ax3, opts, bg, fg, gridC, tSec, LTM_sec, LCM_sec, r_sc, r_earth, iLTM, iLCM, dvHat, unitLabel)
pad = opts.insetDaysPad * 86400;
t1 = min(LTM_sec, LCM_sec) - pad;
t2 = max(LTM_sec, LCM_sec) + pad;
mask = (tSec >= t1) & (tSec <= t2);
idxAll = find(mask);
if numel(idxAll) < 50
    idxAll = max(1,iLTM-200):min(numel(tSec), iLCM+200);
end
idxInset = decimateWithAnchors(numel(tSec), min(opts.maxInsetPoints, numel(idxAll)), [idxAll(1) iLTM iLCM idxAll(end)]);
idxInset = idxInset(idxInset>=idxAll(1) & idxInset<=idxAll(end));
idxPre  = unique([idxInset(idxInset <= iLTM); iLTM]);
idxMid  = unique([idxInset(idxInset >= iLTM & idxInset <= iLCM); iLTM; iLCM]);
idxPost = unique([idxInset(idxInset >= iLCM); iLCM]);
pts = r_sc(idxInset,:);
if opts.insetShowEarth && ~isempty(r_earth)
    pts = [pts; r_earth(idxInset,:)];
end
mins = min(pts,[],1);
maxs = max(pts,[],1);
ctr  = (mins+maxs)/2;
half = max(maxs-mins)/2;
half = max(half, 1e-6);
half = half * opts.insetPad;
pos = ax3.Position;
s   = opts.insetSizeFrac;
w   = pos(3)*s;
h   = pos(4)*s;
m   = 0.02;
switch lower(char(string(opts.insetCorner)))
    case 'northwest'
        x = pos(1) + pos(3)*m;
        y = pos(2) + pos(4) - h - pos(4)*m;
    case 'southeast'
        x = pos(1) + pos(3) - w - pos(3)*m;
        y = pos(2) + pos(4)*m;
    case 'southwest'
        x = pos(1) + pos(3)*m;
        y = pos(2) + pos(4)*m;
    otherwise
        x = pos(1) + pos(3) - w - pos(3)*m;
        y = pos(2) + pos(4) - h - pos(4)*m;
end
axI = axes('Parent',hFig,'Units','normalized','Position',[x y w h]);
hold(axI,'on'); box(axI,'on'); grid(axI,'on');
axis(axI,'equal'); axis(axI,'vis3d');
view(axI, opts.viewVec);
set(axI,'Color',bg,'XColor',fg,'YColor',fg,'ZColor',fg, ...
    'GridColor',gridC,'GridAlpha',0.20, ...
    'FontName',opts.fontName,'FontSize',max(opts.fontSize-1,9), ...
    'LineWidth', 1.0);
if opts.insetShowEarth && ~isempty(r_earth)
    plot3(axI, r_earth(idxInset,1), r_earth(idxInset,2), r_earth(idxInset,3), ':', ...
        'Color', gridC, 'LineWidth', 1.0, 'HandleVisibility','off');
end
plot3(axI, r_sc(idxInset,1), r_sc(idxInset,2), r_sc(idxInset,3), '-', ...
    'Color', gridC, 'LineWidth', 0.9, 'HandleVisibility','off');
plot3(axI, r_sc(idxPre,1), r_sc(idxPre,2), r_sc(idxPre,3), '-', ...
    'Color', opts.cPre, 'LineWidth', 1.3, 'HandleVisibility','off');
plot3(axI, r_sc(idxMid,1), r_sc(idxMid,2), r_sc(idxMid,3), '-', ...
    'Color', opts.cMid, 'LineWidth', 2.6, 'HandleVisibility','off');
plot3(axI, r_sc(idxPost,1), r_sc(idxPost,2), r_sc(idxPost,3), '-', ...
    'Color', opts.cPost, 'LineWidth', 2.2, 'HandleVisibility','off');
plot3(axI, r_sc(iLTM,1), r_sc(iLTM,2), r_sc(iLTM,3), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', fg, 'LineWidth', 1.0);
plot3(axI, r_sc(iLCM,1), r_sc(iLCM,2), r_sc(iLCM,3), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.2 1.0 0.2], 'MarkerEdgeColor', fg, 'LineWidth', 1.0);
boxSizeI = 2*half;
arrowLen = opts.dvFrac * boxSizeI;
r0 = r_sc(iLTM,:);
if norm(r0) > eps, rhat = r0 / norm(r0); else, rhat = dvHat; end
base = r0 + rhat * (opts.dvOffsetFrac * boxSizeI);
quiver3(axI, base(1), base(2), base(3), ...
    dvHat(1)*arrowLen, dvHat(2)*arrowLen, dvHat(3)*arrowLen, 0, ...
    'Color', opts.dvColor, 'LineWidth', 2.0, 'MaxHeadSize', 1.0, 'HandleVisibility','off');
xlim(axI, [ctr(1)-half ctr(1)+half]);
ylim(axI, [ctr(2)-half ctr(2)+half]);
zlim(axI, [ctr(3)-half ctr(3)+half]);
axI.XTickLabel = [];
axI.YTickLabel = [];
axI.ZTickLabel = [];
title(axI, opts.insetTitle, 'Color', fg, 'FontWeight','bold');
text(axI, ctr(1)-half, ctr(2)-half, ctr(3)+half, sprintf('  (%s)', unitLabel), ...
    'Color', fg, 'FontSize', max(opts.fontSize-2,8), 'FontWeight','bold');
end
function optsOut = applyDefaults(optsIn, defIn)
optsOut = optsIn;
fn = fieldnames(defIn);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(optsOut,f) || isempty(optsOut.(f))
        optsOut.(f) = defIn.(f);
    end
end
end
function i = nearestIndex(tVec, t)
[~, i] = min(abs(tVec - t));
i = max(1, min(i, numel(tVec)));
end
function [bg, fg, gridC] = themeColors(theme)
t = lower(char(string(theme)));
switch t
    case 'dark'
        bg = [0.03 0.03 0.05];
        fg = [0.92 0.92 0.95];
        gridC = [0.55 0.55 0.60];
    otherwise
        bg = [1 1 1];
        fg = [0.10 0.10 0.10];
        gridC = [0.70 0.70 0.70];
end
end
function [hFig, ax3, axR] = makeLayout(opts)
hFig = figure('Color',[1 1 1]);
hFig = gcf;
useTiled = (exist('tiledlayout','file') == 2);
if opts.showRangePanel
    if useTiled
        t = tiledlayout(hFig, 2, 1, 'TileSpacing','compact', 'Padding','compact');
        try, t.RowHeight = {'3x','1x'}; catch, end
        ax3 = nexttile(t,1);
        axR = nexttile(t,2);
    else
        ax3 = subplot(2,1,1,'Parent',hFig);
        axR = subplot(2,1,2,'Parent',hFig);
    end
else
    ax3 = axes('Parent', hFig);
    axR = [];
end
end
function idx = decimateWithAnchors(N, maxPts, anchors)
anchors = unique(anchors(:));
anchors = anchors(anchors>=1 & anchors<=N);
if N <= maxPts
    idx = (1:N).';
    return
end
idx = unique(round(linspace(1, N, maxPts))).';
idx = unique([idx; anchors]);
idx = sort(idx);
end
function [xl, yl, zl, boxSize] = computeBoxLimits(pts, opts)
mins = min(pts, [], 1);
maxs = max(pts, [], 1);
span = maxs - mins;
switch lower(char(string(opts.centerMode)))
    case 'sun'
        ctr = [0 0 0];
        half = max(max(abs(pts),[],1));
        half = max(half);
    otherwise
        ctr  = (mins + maxs)/2;
        half = max(span)/2;
end
half = max(half, 1e-6);
half = half * opts.limPad;
xl = [ctr(1)-half, ctr(1)+half];
yl = [ctr(2)-half, ctr(2)+half];
zl = [ctr(3)-half, ctr(3)+half];
boxSize = 2*half;
end
function addStarField(ax, pts, nStars, fg)
mins = min(pts,[],1); maxs = max(pts,[],1);
ctr  = (mins+maxs)/2;
half = max(maxs-mins)/2;
R = 1.3*half;
P = (rand(nStars,3)*2 - 1);
P = P ./ vecnorm(P,2,2);
P = ctr + P * R;
h = scatter3(ax, P(:,1),P(:,2),P(:,3), 6, fg, 'filled', 'HandleVisibility','off');
try
    h.MarkerFaceAlpha = 0.12;
    h.MarkerEdgeAlpha = 0.0;
catch
end
end