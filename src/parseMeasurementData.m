function meas = parseMeasurementData(filename, t0_et)
    if nargin < 2
        t0_et = 0;
    end
    T = readtable(filename, 'VariableNamingRule','preserve');
    varNames = string(T.Properties.VariableNames);
    timeVar = "";
    candAbs = ["Time (sec past J2000 - ET)", "Time_ET", "Time_ET_s"];
    candRel = ["Time (s)", "Time_rel_s", "Time since detection (s)"];
    for c = candAbs
        if any(varNames == c)
            timeVar = c;
            isAbsET = true;
            break;
        end
    end
    if timeVar == ""
        for c = candRel
            if any(varNames == c)
                timeVar = c;
                isAbsET = false;
                break;
            end
        end
    end
    if timeVar == ""
        error('parseMeasurementData:NoTimeColumn', ...
              'No recognizable time column found in %s', filename);
    end
    t_raw = T.(timeVar);
    if isAbsET
        if t0_et ~= 0
            warning('parseMeasurementData:IgnoringT0', ...
                'Time column appears to be absolute ET; t0_et is ignored.');
        end
        t_et = t_raw;
    else
        if t0_et == 0
            warning('parseMeasurementData:MissingT0', ...
                ['Time column appears relative to detection but t0_et=0. ' ...
                 'Make sure this is what you intend.']);
        end
        t_et = t0_et + t_raw;
    end
    if any(varNames == "Station ID")
        station = T.("Station ID");
    else
        error('parseMeasurementData:NoStationID', ...
              'No "Station ID" column found in %s', filename);
    end
    rrVar = "";
    candRR = ["Range Rate (km/s)", "Range rate (km/s)", "RangeRate_km_s"];
    for c = candRR
        if any(varNames == c)
            rrVar = c;
            break;
        end
    end
    if rrVar == ""
        error('parseMeasurementData:NoRangeRate', ...
              'No range-rate column found in %s', filename);
    end
    rr = T.(rrVar);
    rangeVar = "";
    candR = ["Range (km)", "Range_km"];
    for c = candR
        if any(varNames == c)
            rangeVar = c;
            break;
        end
    end
    if rangeVar ~= ""
        range = T.(rangeVar);
    else
        range = nan(size(rr));
    end
    [t_sorted, idx] = sort(t_et);
    meas.t         = t_sorted(:);
    meas.stationID = station(idx);
    meas.range     = range(idx);
    meas.rr        = rr(idx);
end