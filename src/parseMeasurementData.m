function meas = parseMeasurementData(filename, t0_et)
% parseMeasurementData  Read truth measurement CSV and return ET timestamps.
%
% Inputs:
%   filename - CSV file name
%   t0_et    - (optional) detection epoch in ET [sec past J2000]
%              Used only if the file time is relative to detection.
%
% Output:
%   meas struct with fields:
%       .t         - Nx1 measurement times in ET [s]
%       .stationID - Nx1 station IDs (1..4)
%       .range     - Nx1 range [km]   (NaN if not present in file)
%       .rr        - Nx1 range-rate [km/s]

    if nargin < 2
        t0_et = 0;  % default: assume times are already ET
    end

    T = readtable(filename, 'VariableNamingRule','preserve');
    varNames = string(T.Properties.VariableNames);

    % -------- Time column (ET vs relative) -------- %
    timeVar = "";
    % Candidates for absolute ET
    candAbs = ["Time (sec past J2000 - ET)", "Time_ET", "Time_ET_s"];
    % Candidates for relative-to-detection
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
        % File already in ET seconds past J2000
        if t0_et ~= 0
            warning('parseMeasurementData:IgnoringT0', ...
                'Time column appears to be absolute ET; t0_et is ignored.');
        end
        t_et = t_raw;
    else
        % File is relative to detection; need t0_et
        if t0_et == 0
            warning('parseMeasurementData:MissingT0', ...
                ['Time column appears relative to detection but t0_et=0. ' ...
                 'Make sure this is what you intend.']);
        end
        t_et = t0_et + t_raw;
    end

    % -------- Station ID -------- %
    if any(varNames == "Station ID")
        station = T.("Station ID");
    else
        error('parseMeasurementData:NoStationID', ...
              'No "Station ID" column found in %s', filename);
    end

    % -------- Range-rate column -------- %
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

    % -------- Optional range column -------- %
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
        range = nan(size(rr));  % no range in 0–6D file
    end

    % -------- Pack into struct and sort by time -------- %
    [t_sorted, idx] = sort(t_et);

    meas.t         = t_sorted(:);
    meas.stationID = station(idx);
    meas.range     = range(idx);
    meas.rr        = rr(idx);
end
