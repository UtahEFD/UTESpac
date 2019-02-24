function output = avg(data,info,tableNames,output,headers,sensorInfo)

for i = 1:length(data)
    display(sprintf('\naveraging %s',tableNames{i}))
    currentTable = data{i};
    if ~isempty(currentTable)
        
        % find table header
        header = headers{i};
        
        % find table frequency
        tableFreq = 1/((currentTable(2,1) - currentTable(1,1))*(3600*24));
        
        % find number of samples per averaging period
        numSamplesPerPeriod = round(info.avgPer*60*tableFreq);
        if numSamplesPerPeriod<1
            numSamplesPerPeriod = 1;
            warning(['Averaging Period of ', num2str(info.avgPer),...
                ' is too short for table: ', tableNames{i}]);
        end
        
        % find number of averaging periods
        numAvgPeriods = size(currentTable,1)/numSamplesPerPeriod;
        
        
        % preallocate avgTable
        avgTable = nan(numAvgPeriods,size(currentTable,2));
        
        % iterate through avgTable rows
        for j = 1:size(avgTable,1)
            
            % store local data
            localData = currentTable((j-1)*numSamplesPerPeriod+1:j*numSamplesPerPeriod,:);
            
            % average current data and store in avgTable
            avgTable(j,:) = nanmean(localData);
            
            % date is last value of averaging period
            avgTable(j,1) = localData(end,1);  
            
            % check for wind direction from wind bird and find mean unit vector direction
            if isfield(sensorInfo,'birdDir')
                for k = 1:length(sensorInfo.birdDir(:,1))
                    if sensorInfo.birdDir(k,1) == i
                        WS = localData(:,sensorInfo.birdSpd(k,2));
                        WD = localData(:,sensorInfo.birdDir(k,2)).*pi./180;
                        v = nanmean(WS.*sin(WD));
                        u = nanmean(WS.*cos(WD));
                        avgTable(j,sensorInfo.birdDir(k,2)) = mod(atan2(v,u)*180/pi,360);
                    end
                end
            end
            
        end
        tableName = tableNames{i};
        output.(tableName) = avgTable;
        output.([tableName,'Header']) = header;
    end
end
end