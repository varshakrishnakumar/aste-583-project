function writeMissionFollowCamVideo(T, r_sc, r_body2, ET0, LTM, LCM, deltaVLTM, opts)
if nargin < 8 || isempty(opts), opts = struct(); end
def.videoFile      = 'mission_followcam.mp4';
def.fps            = 30;
def.videoSeconds   = 60;
def.nFrames         = [];
def.figSize        = [100 100 1280 720];
def.theme          = 'dark';
def.showRangePanel = true;
def.showFullPath       = true;
def.maxFullPathPoints  = 12000;
def.tailPoints         = 900;
def.followMode     = 'axis';
def.viewVec        = [-135 25];
def.camViewAngle   = 12;
def.camBackFactor  = 2.2;
def.camUpFactor    = 0.8;
def.camSmoothing   = 0.12;
def.followRangeKm      = [];
def.followRangeFactor  = 6;
def.minRangeKm         = 2e4;
def.maxRangeKm         = 4e6;
def.primaryName    = 'Earth';
def.body2Name      = 'Moon';
def.body3Name      = 'Moon';
def.showBody3      = false;
def.r_body3        = [];
def.msBody         = 70;
def.msSC           = 55;
def.cPrimary       = [0.2 0.6 1.0];
def.cBody2         = [0.45 0.45 0.45];
def.cBody3         = [0.70 0.70 0.70];
def.cPre      = [0.65 0.65 0.65];
def.cLTM_LCM  = [0.85 0.33 0.10];
def.cLCM_CA   = [0.00 0.45 0.74];
def.cPost     = [0.30 0.75 0.35];
def.lwTail    = 2.4;
def.lwFaint   = 1.2;
def.dvFrac        = 0.18;
def.dvColor       = [0.55 0.00 0.65];
def.dvLW          = 2.4;
def.dvHead        = 1.4;
def.dvOffsetFrac  = 0.12;
def.dvShowSeconds = 2*3600;
def.v_sc = [];
opts = applyDefaults(opts, def);
T = T(:);
r_sc = r_sc(:,1:3);
N = size(r_sc,1);
r_body2 = normalizeBodyTraj(r_body2, N);
if opts.showBody3 && ~isempty(opts.r_body3)
    r_body3 = normalizeBodyTraj(opts.r_body3, N);
else
    r_body3 = [];
end
iLTM = nearestIndex(T, LTM);
iLCM = nearestIndex(T, LCM);
if ~isempty(r_body3)
    dCA = vecnorm(r_sc - r_body3, 2, 2);
else
    dCA = vecnorm(r_sc - r_body2, 2, 2);
end
[dmin, iCA] = min(dCA);
tDays = (T - ET0)/86400;
rPrim = vecnorm(r_sc, 2, 2);
r2    = vecnorm(r_sc - r_body2, 2, 2);
if ~isempty(r_body3)
    r3 = vecnorm(r_sc - r_body3, 2, 2);
else
    r3 = [];
end
dv = deltaVLTM(:);
dvMag_mps = 1e3 * norm(dv);
if ~isempty(opts.nFrames)
    nFrames = max(2, round(opts.nFrames));
else
    nFrames = max(2, round(opts.fps * opts.videoSeconds));
end
idxFrame = unique(round(linspace(1, N, nFrames)));
nFrames  = numel(idxFrame);
[hFig, ax3, axR] = makeLayout(opts);
set(hFig, 'Color', [1 1 1], 'Renderer', 'opengl', 'InvertHardcopy', 'off');
set(hFig, 'Position', opts.figSize);
[bg, fg, gridC] = themeColors(opts.theme);
set(hFig, 'Color', bg);
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
axis(ax3,'equal'); axis(ax3,'vis3d');
view(ax3, opts.viewVec);
set(ax3, 'Color', bg, ...
    'XColor', fg, 'YColor', fg, 'ZColor', fg, ...
    'GridColor', gridC, 'GridAlpha', 0.25);
xlabel(ax3, sprintf('X_{%s} (km)', opts.primaryName));
ylabel(ax3, sprintf('Y_{%s} (km)', opts.primaryName));
zlabel(ax3, sprintf('Z_{%s} (km)', opts.primaryName));
half0 = computeHalfRange(1, r2, r3, opts);
setAxisWindow(ax3, r_sc(1,:), half0);
scatter3(ax3, 0,0,0, opts.msBody, 'filled', ...
    'MarkerFaceColor', opts.cPrimary, 'MarkerEdgeColor', fg, ...
    'DisplayName', opts.primaryName);
text(ax3, 0,0,0, ['  ' opts.primaryName], 'Color', fg, 'FontWeight','bold');
hB2 = scatter3(ax3, r_body2(1,1), r_body2(1,2), r_body2(1,3), opts.msBody, 'filled', ...
    'MarkerFaceColor', opts.cBody2, 'MarkerEdgeColor', fg, ...
    'DisplayName', opts.body2Name);
hLineB2 = plot3(ax3, [0 r_body2(1,1)], [0 r_body2(1,2)], [0 r_body2(1,3)], ':', ...
    'Color', gridC, 'LineWidth', 1.1, 'HandleVisibility','off');
text(ax3, r_body2(1,1), r_body2(1,2), r_body2(1,3), ['  ' opts.body2Name], 'Color', fg, 'FontWeight','bold');
if ~isempty(r_body3)
    hB3 = scatter3(ax3, r_body3(1,1), r_body3(1,2), r_body3(1,3), opts.msBody, 'filled', ...
        'MarkerFaceColor', opts.cBody3, 'MarkerEdgeColor', fg, ...
        'DisplayName', opts.body3Name);
    hLineB3 = plot3(ax3, [0 r_body3(1,1)], [0 r_body3(1,2)], [0 r_body3(1,3)], ':', ...
        'Color', gridC, 'LineWidth', 1.0, 'HandleVisibility','off');
    text(ax3, r_body3(1,1), r_body3(1,2), r_body3(1,3), ['  ' opts.body3Name], 'Color', fg, 'FontWeight','bold');
else
    hB3 = []; hLineB3 = [];
end
if opts.showFullPath
    idxFull = unique(round(linspace(1, N, min(N, opts.maxFullPathPoints))));
    cFull = mixColor(fg, bg, 0.65);
    plot3(ax3, r_sc(idxFull,1), r_sc(idxFull,2), r_sc(idxFull,3), '-', ...
        'Color', cFull, 'LineWidth', opts.lwFaint, 'DisplayName', 'Trajectory (full)');
end
hPre  = plot3(ax3, NaN,NaN,NaN, '-', 'Color', opts.cPre,     'LineWidth', opts.lwFaint, 'HandleVisibility','off');
h12   = plot3(ax3, NaN,NaN,NaN, '-', 'Color', opts.cLTM_LCM, 'LineWidth', opts.lwTail,  'HandleVisibility','off');
h23   = plot3(ax3, NaN,NaN,NaN, '-', 'Color', opts.cLCM_CA,  'LineWidth', opts.lwTail,  'HandleVisibility','off');
hPost = plot3(ax3, NaN,NaN,NaN, '-', 'Color', opts.cPost,    'LineWidth', opts.lwFaint, 'HandleVisibility','off');
hSC = plot3(ax3, r_sc(1,1), r_sc(1,2), r_sc(1,3), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', fg, ...
    'LineWidth', 1.2, 'DisplayName', 'Spacecraft');
plot3(ax3, r_sc(iLTM,1), r_sc(iLTM,2), r_sc(iLTM,3), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', fg, 'LineWidth', 1.2, ...
    'DisplayName','LTM');
plot3(ax3, r_sc(iLCM,1), r_sc(iLCM,2), r_sc(iLCM,3), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 1.0 0.2], 'MarkerEdgeColor', fg, 'LineWidth', 1.2, ...
    'DisplayName','LCM');
plot3(ax3, r_sc(iCA,1), r_sc(iCA,2), r_sc(iCA,3), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.4, ...
    'DisplayName','CA');
dvHat = dv(:).';
if any(~isfinite(dvHat)) || norm(dvHat) < eps
    dvHat = [1 0 0];
else
    dvHat = dvHat / norm(dvHat);
end
halfLTM = computeHalfRange(iLTM, r2, r3, opts);
arrowLen = opts.dvFrac * (2*halfLTM);
r0 = r_sc(iLTM,:);
if norm(r0) > eps
    rhat = r0 / norm(r0);
else
    rhat = dvHat;
end
base = r0 + rhat * (opts.dvOffsetFrac * (2*halfLTM));
hDV = quiver3(ax3, base(1), base(2), base(3), ...
    dvHat(1)*arrowLen, dvHat(2)*arrowLen, dvHat(3)*arrowLen, 0, ...
    'Color', opts.dvColor, 'LineWidth', opts.dvLW, 'MaxHeadSize', opts.dvHead, ...
    'DisplayName','\DeltaV @ LTM', 'Visible','off');
title(ax3, sprintf('Follow-cam mission view | CA min = %.1f km', dmin), 'Color', fg, 'FontWeight','bold');
legend(ax3, 'Location','northeastoutside', 'TextColor', fg);
hBox = annotation(hFig, 'textbox', [0.02 0.72 0.42 0.24], 'String', '', ...
    'FitBoxToText','on', 'BackgroundColor', bg, 'EdgeColor', gridC, ...
    'Color', fg, 'FontName', 'Helvetica', 'FontSize', 10);
hNow = [];
if opts.showRangePanel && ~isempty(axR) && isvalid(axR)
    hold(axR,'on'); grid(axR,'on'); box(axR,'on');
    set(axR, 'Color', bg, 'XColor', fg, 'YColor', fg, 'GridColor', gridC, 'GridAlpha', 0.25);
    semilogy(axR, tDays, rPrim, '-',  'Color', mixColor(fg,bg,0.35), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('|r_{SC-%s}|', opts.primaryName));
    semilogy(axR, tDays, r2,   '--', 'Color', opts.cLCM_CA, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('|r_{SC-%s}|', opts.body2Name));
    if ~isempty(r3)
        semilogy(axR, tDays, r3, ':', 'Color', opts.cBody3, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('|r_{SC-%s}|', opts.body3Name));
    end
    addXLine(axR, tDays(iLTM), ':', 'LTM', [1 0.2 0.2]);
    addXLine(axR, tDays(iLCM), ':', 'LCM', [0.2 0.8 0.2]);
    addXLine(axR, tDays(iCA),  ':', 'CA',  fg);
    xlabel(axR, 'Days since detection');
    ylabel(axR, 'Range (km) (log)');
    legend(axR,'Location','northeast');
    hNow = addXLine(axR, tDays(1), '-', 'Now', fg);
end
vw = VideoWriter(opts.videoFile, 'MPEG-4');
vw.FrameRate = opts.fps;
try, vw.Quality = 95; catch, end
open(vw);
ctrSm  = r_sc(idxFrame(1),:);
halfSm = half0;
for k = 1:nFrames
    i = idxFrame(k);
    set(hB2, 'XData', r_body2(i,1), 'YData', r_body2(i,2), 'ZData', r_body2(i,3));
    set(hLineB2, 'XData', [0 r_body2(i,1)], 'YData', [0 r_body2(i,2)], 'ZData', [0 r_body2(i,3)]);
    if ~isempty(r_body3)
        set(hB3, 'XData', r_body3(i,1), 'YData', r_body3(i,2), 'ZData', r_body3(i,3));
        set(hLineB3, 'XData', [0 r_body3(i,1)], 'YData', [0 r_body3(i,2)], 'ZData', [0 r_body3(i,3)]);
    end
    i0 = max(1, i - opts.tailPoints + 1);
    win = i0:i;
    idxPre  = win(win <= iLTM);
    idx12   = win(win >= iLTM & win <= iLCM);
    idx23   = win(win >= iLCM & win <= iCA);
    idxPost = win(win >= iCA);
    setLineXYZ(hPre,  r_sc, idxPre);
    setLineXYZ(h12,   r_sc, idx12);
    setLineXYZ(h23,   r_sc, idx23);
    setLineXYZ(hPost, r_sc, idxPost);
    set(hSC, 'XData', r_sc(i,1), 'YData', r_sc(i,2), 'ZData', r_sc(i,3));
    half = computeHalfRange(i, r2, r3, opts);
    ctr = r_sc(i,:);
    if opts.camSmoothing > 0
        a = max(0, min(1, opts.camSmoothing));
        ctrSm  = (1-a)*ctrSm  + a*ctr;
        halfSm = (1-a)*halfSm + a*half;
    else
        ctrSm  = ctr;
        halfSm = half;
    end
    setAxisWindow(ax3, ctrSm, halfSm);
    if strcmpi(string(opts.followMode), "chase")
        vhat = estimateVelHat(i, T, r_sc, opts.v_sc);
        back = opts.camBackFactor * halfSm;
        up   = opts.camUpFactor   * halfSm;
        campos(ax3, ctrSm - vhat*back + [0 0 1]*up);
        camtarget(ax3, ctrSm);
        camup(ax3, [0 0 1]);
        camva(ax3, opts.camViewAngle);
    end
    if abs(T(i) - LTM) <= opts.dvShowSeconds/2
        set(hDV, 'Visible', 'on');
    else
        set(hDV, 'Visible', 'off');
    end
    phase = phaseName(i, iLTM, iLCM, iCA);
    if isempty(r3)
        s = sprintf([ ...
            't = %.2f d (%s)\n' ...
            '|r_{SC-%s}| = %.3f Mm\n' ...
            '|r_{SC-%s}| = %.3f Mm\n' ...
            'LTM @ %.2f d | LCM @ %.2f d | CA @ %.2f d\n' ...
            '\\DeltaV(LTM) = %.0f m/s'], ...
            tDays(i), phase, ...
            opts.primaryName, rPrim(i)/1e6, ...
            opts.body2Name,   r2(i)/1e6, ...
            tDays(iLTM), tDays(iLCM), tDays(iCA), ...
            dvMag_mps);
    else
        s = sprintf([ ...
            't = %.2f d (%s)\n' ...
            '|r_{SC-%s}| = %.3f Mm\n' ...
            '|r_{SC-%s}| = %.3f Mm\n' ...
            '|r_{SC-%s}| = %.3f Mm\n' ...
            'LTM @ %.2f d | LCM @ %.2f d | CA @ %.2f d\n' ...
            '\\DeltaV(LTM) = %.0f m/s'], ...
            tDays(i), phase, ...
            opts.primaryName, rPrim(i)/1e6, ...
            opts.body2Name,   r2(i)/1e6, ...
            opts.body3Name,   r3(i)/1e6, ...
            tDays(iLTM), tDays(iLCM), tDays(iCA), ...
            dvMag_mps);
    end
    set(hBox, 'String', s);
    if ~isempty(hNow)
        updateXLine(hNow, tDays(i), axR);
    end
    drawnow limitrate nocallbacks;
    writeVideo(vw, getframe(hFig));
end
close(vw);
disp(['Wrote video: ' char(opts.videoFile)]);
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
function i = nearestIndex(Tvec, t)
[~, i] = min(abs(Tvec - t));
i = max(1, min(i, numel(Tvec)));
end
function [bg, fg, gridC] = themeColors(theme)
t = lower(char(string(theme)));
switch t
    case 'dark'
        bg = [0 0 0];
        fg = [0.92 0.92 0.92];
        gridC = [0.55 0.55 0.55];
    otherwise
        bg = [1 1 1];
        fg = [0.10 0.10 0.10];
        gridC = [0.70 0.70 0.70];
end
end
function c = mixColor(c1, c2, alpha)
alpha = max(0, min(1, alpha));
c = (1-alpha)*c1 + alpha*c2;
end
function rN = normalizeBodyTraj(r, N)
r = r(:,1:3);
if size(r,1) == 1
    rN = repmat(r, N, 1);
elseif size(r,1) ~= N
    error('Body trajectory must be Nx3 (same N as r_sc) or 1x3 constant.');
else
    rN = r;
end
end
function [hFig, ax3, axR] = makeLayout(opts)
hFig = figure('Color',[1 1 1]);
useTiled = (exist('tiledlayout','file') == 2);
if opts.showRangePanel
    if useTiled
        t = tiledlayout(hFig, 2, 1, ...
            'TileSpacing','compact', ...
            'Padding','compact');
        ax3 = nexttile(t,1);
        axR = nexttile(t,2);
        try
            t.RowHeight = {'1x','1x'};
        catch
            ax3.Position = [0.08 0.55 0.85 0.40];
            axR.Position = [0.08 0.08 0.85 0.40];
        end
    else
        ax3 = subplot(2,1,1,'Parent',hFig);
        axR = subplot(2,1,2,'Parent',hFig);
        ax3.Position = [0.08 0.55 0.85 0.40];
        axR.Position = [0.08 0.08 0.85 0.40];
    end
else
    ax3 = axes('Parent', hFig);
    axR = [];
end
end
function half = computeHalfRange(i, r2, r3, opts)
if ~isempty(opts.followRangeKm)
    half = opts.followRangeKm;
else
    if isempty(r3)
        localScale = r2(i);
    else
        localScale = min(r2(i), r3(i));
    end
    half = opts.followRangeFactor * localScale;
    half = max(opts.minRangeKm, min(opts.maxRangeKm, half));
end
end
function setAxisWindow(ax, ctr, half)
xlim(ax, ctr(1) + [-half half]);
ylim(ax, ctr(2) + [-half half]);
zlim(ax, ctr(3) + [-half half]);
end
function setLineXYZ(h, r, idx)
if isempty(idx)
    set(h,'XData',NaN,'YData',NaN,'ZData',NaN);
else
    set(h,'XData',r(idx,1),'YData',r(idx,2),'ZData',r(idx,3));
end
end
function vhat = estimateVelHat(i, T, r, v_sc)
N = size(r,1);
if ~isempty(v_sc)
    v = v_sc(i,:);
else
    if i <= 1
        dr = r(2,:) - r(1,:); dt = T(2)-T(1);
    elseif i >= N
        dr = r(N,:) - r(N-1,:); dt = T(N)-T(N-1);
    else
        dr = r(i+1,:) - r(i-1,:); dt = T(i+1)-T(i-1);
    end
    dt = max(dt, eps);
    v = dr / dt;
end
nv = norm(v);
if nv < eps
    vhat = [1 0 0];
else
    vhat = v / nv;
end
end
function txt = phaseName(i, iLTM, iLCM, iCA)
if i < iLTM
    txt = 'Detect \rightarrow LTM';
elseif i < iLCM
    txt = 'LTM \rightarrow LCM';
elseif i < iCA
    txt = 'LCM \rightarrow CA';
else
    txt = 'Post-CA';
end
end
function h = addXLine(ax, x, ls, labelStr, c)
if exist('xline','file') == 2
    h = xline(ax, x, ls, labelStr, 'Color', c, 'LineWidth', 1.3, 'LabelOrientation','horizontal');
else
    yl = ylim(ax);
    h = plot(ax, [x x], yl, ls, 'Color', c, 'LineWidth', 1.3);
    if ~isempty(labelStr)
        text(ax, x, yl(2), [' ' labelStr], 'Color', c, 'VerticalAlignment','top');
    end
end
end
function updateXLine(h, x, ax)
    if isprop(h, 'Value')
        h.Value = x;
    else
        yl = ylim(ax);
        set(h, 'XData', [x x], 'YData', yl);
    end
end