function measData = parseMeasurementData(filename)
% parseMeasurementData Reads Lunar Trailblazer measurement CSV file
%   Input:
%     filename - string path to the CSV file (e.g., 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv')
%   Output:
%     measData - struct with parsed fields: time (ET), stationID, range, rangeRate

    % Read CSV assuming it has a header row
    raw = readtable(filename);

    % Rename fields for consistency
    measData.time        = raw.time_et;        % seconds past J2000
    measData.stationID   = raw.station_id;     % station number (1–4)
    measData.range       = raw.range_km;       % km
    measData.rangeRate   = raw.range_rate_kms; % km/s
end
