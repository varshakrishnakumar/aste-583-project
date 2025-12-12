function hFig = plotMissionTrajectory3D(Ttotal, XrelEarth, xMoon, LTM, LCM, deltaVLTM, opts)
if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts,'bodyScale'),      opts.bodyScale      = 200; end
if ~isfield(opts,'earthRadius_km'), opts.earthRadius_km = 6371; end
if ~isfield(opts,'moonRadius_km'),  opts.moonRadius_km  = 1737; end
if ~isfield(opts,'earthTex'),       opts.earthTex       = '';  end
if ~isfield(opts,'moonTex'),        opts.moonTex        = '';  end
if ~isfield(opts,'arrowFrac'),      opts.arrowFrac      = 0.02; end
if ~isfield(opts,'viewVec'),        opts.viewVec        = [-135 25]; end
T = Ttotal(:);
r_sc   = XrelEarth(:,1:3);
r_moon = xMoon(:,1:3);
[~, iLTM] = min(abs(T - LTM));
[~, iLCM] = min(abs(T - LCM));
spanXYZ = max(r_sc,[],1) - min(r_sc,[],1);
span    = max(spanXYZ);
if span <= 0, span = 1; end
arrowLen = opts.arrowFrac * span;
Re = opts.earthRadius_km * opts.bodyScale;
Rm = opts.moonRadius_km  * opts.bodyScale;
hFig = figure('Color','w');
ax = axes('Parent',hFig); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'X_E (km)'); ylabel(ax,'Y_E (km)'); zlabel(ax,'Z_E (km)');
title(ax,'CONOPS 3D: Earth-centered trajectory');
[XS,YS,ZS] = sphere(50);
XS = Re*XS; YS = Re*YS; ZS = Re*ZS;
hE = surf(ax, XS,YS,ZS, 'EdgeColor','none');
if ~isempty(opts.earthTex) && exist(opts.earthTex,'file')
    tex = flipud(imread(opts.earthTex));
    set(hE,'FaceColor','texturemap','CData',tex);
else
    set(hE,'FaceColor',[0.2 0.4 0.9]);
end
plot3(ax, r_moon(:,1), r_moon(:,2), r_moon(:,3), ...
      '--', 'LineWidth', 1.2, 'DisplayName','Moon path');
[XM,YM,ZM] = sphere(40);
XM = Rm*XM + r_moon(iLCM,1);
YM = Rm*YM + r_moon(iLCM,2);
ZM = Rm*ZM + r_moon(iLCM,3);
hM = surf(ax, XM,YM,ZM, 'EdgeColor','none');
if ~isempty(opts.moonTex) && exist(opts.moonTex,'file')
    texM = flipud(imread(opts.moonTex));
    set(hM,'FaceColor','texturemap','CData',texM);
else
    set(hM,'FaceColor',[0.7 0.7 0.7]);
end
hTraj = plot3(ax, r_sc(:,1), r_sc(:,2), r_sc(:,3), ...
              'b-', 'LineWidth', 1.5, 'DisplayName','Spacecraft');
hLTM = plot3(ax, r_sc(iLTM,1), r_sc(iLTM,2), r_sc(iLTM,3), ...
             'ro','MarkerSize',8,'LineWidth',1.5,'DisplayName','LTM');
hLCM = plot3(ax, r_sc(iLCM,1), r_sc(iLCM,2), r_sc(iLCM,3), ...
             'go','MarkerSize',8,'LineWidth',1.5,'DisplayName','LCM');
dv = deltaVLTM(:).';
if norm(dv) > 0
    dv_hat = dv / norm(dv);
else
    dv_hat = [1 0 0];
end
hDV = quiver3(ax, r_sc(iLTM,1), r_sc(iLTM,2), r_sc(iLTM,3), ...
              dv_hat(1)*arrowLen, dv_hat(2)*arrowLen, dv_hat(3)*arrowLen, 0, ...
              'Color','r','LineWidth',2,'MaxHeadSize',0.6,'DisplayName','\DeltaV at LTM');
pad = 0.03*span;
xlim(ax,[min(r_sc(:,1))-pad, max(r_sc(:,1))+pad]);
ylim(ax,[min(r_sc(:,2))-pad, max(r_sc(:,2))+pad]);
zlim(ax,[min(r_sc(:,3))-pad, max(r_sc(:,3))+pad]);
view(ax, opts.viewVec);
lighting(ax,'gouraud'); material(ax,'dull');
camlight(ax,'headlight');
legend(ax,[hTraj hLTM hLCM hDV],'Location','northeast');
end