function makeResidualPlots(t_meas_hr, prefit, postfit, Niter, maxPerFig, baseTitle, yLabelUnit, prefix)
    
    figIndex = 1;
    tileIndex = 1;

    figure;
    tlo = tiledlayout(maxPerFig,1);

    for it = 1:Niter
        
        % --- Prefit ---
        nexttile; hold on; grid on;
        scatter(t_meas_hr, prefit(it,:), 10, 'filled');
        xline(24, 'k--', 'DCO', 'LabelVerticalAlignment','bottom');
        ylabel(sprintf('Iter %d', it));
        title(sprintf('Prefit %s residuals [%s]', baseTitle, yLabelUnit));
        tileIndex = tileIndex + 1;

        if tileIndex > maxPerFig
            fname = sprintf('%s_prefit_postfit_fig_%02d.png', prefix, figIndex);
            saveas(gcf, fname);
            fprintf('Saved: %s\n', fname);

            figIndex = figIndex + 1;
            figure; tlo = tiledlayout(maxPerFig,1);
            tileIndex = 1;
        end

        % --- Postfit ---
        nexttile; hold on; grid on;
        scatter(t_meas_hr, postfit(it,:), 10, 'filled');
        xline(24, 'k--', 'DCO', 'LabelVerticalAlignment','bottom');
        title(sprintf('Postfit %s residuals [%s]', baseTitle, yLabelUnit));
        tileIndex = tileIndex + 1;

        if tileIndex > maxPerFig
            fname = sprintf('%s_prefit_postfit_fig_%02d.png', prefix, figIndex);
            saveas(gcf, fname);
            fprintf('Saved: %s\n', fname);

            figIndex = figIndex + 1;
            figure; tlo = tiledlayout(maxPerFig,1);
            tileIndex = 1;
        end
    end

    fname = sprintf('%s_prefit_postfit_fig_%02d.png', prefix, figIndex);
    saveas(gcf, fname);
    fprintf('Saved: %s\n', fname);
end
