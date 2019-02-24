function output = windStats(output,sensorInfo,tableNames,info)

% find number sonic info
if isfield(sensorInfo,'u')
    numSonics = size(sensorInfo.u,1);
    towerBearing = info.tower;
    windEnvelope = info.windDirectionTest.envelopeSize;
    
    % create warning field within output to store warning messages
    output.warnings = cell(1);
    
    for i = 1:numSonics
        try
            tble = sensorInfo.u(i,1);
            bearing = sensorInfo.u(i,4);
            sonHeight = sensorInfo.u(i,3);
            uCol = sensorInfo.u(sensorInfo.u(:,3)==sonHeight,2);
            vCol = sensorInfo.v(sensorInfo.v(:,3)==sonHeight,2);
            tableName = tableNames{tble};
            % check sonic manufacturere
            if sensorInfo.u(i,5)
                u = output.(tableName)(:,uCol);
                v = output.(tableName)(:,vCol);
            else % for RMYoung swap u and v and multiply v by -1
                u = output.(tableName)(:,vCol);
                v = output.(tableName)(:,uCol)*-1;
            end
            t = output.(tableName)(:,1);
            
            % find direction and speed
            dir = mod(atan2(-1*v,u)*180/pi+bearing,360);
            spd = (u.^2+v.^2).^0.5;
            flag = zeros(size(dir));
            
            % flag mean wind speed that falls within envelope
            if towerBearing+windEnvelope > 360
                minAngle = towerBearing - windEnvelope;
                maxAngle = towerBearing-360+windEnvelope;
                flag(dir>minAngle) = 1;
                flag(dir<maxAngle) = 1;
            elseif towerBearing-windEnvelope < 0
                minAngle = 360+towerBearing - windEnvelope;
                maxAngle = towerBearing+windEnvelope;
                flag(dir>minAngle) = 1;
                flag(dir<maxAngle) = 1;
            else
                minAngle = towerBearing - windEnvelope;
                maxAngle = towerBearing + windEnvelope;
                temp1 = zeros(size(dir));
                temp2 = zeros(size(dir));
                temp1(dir>minAngle) = 1;
                temp2(dir<maxAngle) = 1;
                flag((temp1+temp2) == 2) = 1;
            end
            output.spdAndDir(1:length(spd),1) = t;
            output.spdAndDir(1:length(spd),i*3-1:i*3+1) = [dir spd flag];
            output.spdAndDirHeader{1} = 'timeStamp';
            output.spdAndDirHeader{i*3-1} = sprintf('%gm direction',sonHeight);
            output.spdAndDirHeader{i*3} = sprintf('%gm speed',sonHeight);
            output.spdAndDirHeader{i*3+1} = sprintf('%gm flag %g<dir<%g',sonHeight,minAngle,maxAngle);
        catch err
            message = strcat(err.message,'@ line',num2str(err.stack.line),' Unable to find wind direction and speed at:',num2str(sonHeight));
            warning(message)
            if isempty(output.warnings{1})
                output.warnings{1,1} = message;
            else output.warnings{end+1,1} = message;
            end
        end
    end
end
