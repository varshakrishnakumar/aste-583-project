function meas = parseMeasurementData(filename, t0_et)
% parseMeasurementData  Read truth measurement CSV and add ET timestamps.
%
% Inputs:
%   filename - CSV file name (e.g. 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv')
%   t0_et    - detection epoch in ET (sec past J2000)
%
% Output:
%   meas struct with fields:
%       .t      - Nx1 measurement times in ET [s]
%       .stationID - Nx1 station IDs (1..4)
%       .range     - Nx1 range [km]   (NaN if not present in file)
%       .rr        - Nx1 range-rate [km/s]

    T = readtable(filename, 'VariableNamingRule','preserve');

    % Required fields
    t_rel    = T.("Time (s)");       % seconds since detection
    station  = T.("Station ID");

    % Range-rate column (handle small naming variations)
    if any(strcmp(T.Properties.VariableNames, "Range Rate (km/s)"))
        rr = T.("Range Rate (km/s)");
    elseif any(strcmp(T.Properties.VariableNames, "Range rate (km/s)"))
        rr = T.("Range rate (km/s)");
    else
        error('No range-rate column found in %s', filename);
    end

    % Optional range column
    if any(strcmp(T.Properties.VariableNames, "Range (km)"))
        range = T.("Range (km)");
    else
        range = nan(size(rr));  % no range in 0-6D file
    end

    meas.t         = t0_et + t_rel;  % convert to ET
    meas.stationID = station;
    meas.range     = range;
    meas.rr        = rr;
end
