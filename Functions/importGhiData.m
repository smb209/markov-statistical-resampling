function GHITable = importGhiData(csv_file_path)
    %% Imports a csv file with the GHI data and returns a timetable
    opts = detectImportOptions(csv_file_path);
    D = readtable(csv_file_path,opts); % <-- | date | time | GHI
    
    % Select the data columns
    dates = D{:,1};
    times = D{:,2};
    GHI   = D{:,3};
    
    %NOTE: there was an extra (redundant time) at 5/17/2017, 7:32 AM.
    indices = find(GHI(:)<0);
    GHI(indices)=0;
    
    % There are some spurious mesurements, saturate this at an average negative
    % value
    GHI(GHI < -15) = -15;
    
    % Convert PV data to a timetable
    tmpTimes = datetime(times,'InputFormat','HH:mm');
    dateTimeVec = dates + hours(hour(tmpTimes)) + minutes(minute(tmpTimes)) + years(2000);
    GHITable = timetable(dateTimeVec,GHI);
end