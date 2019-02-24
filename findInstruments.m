function [sensorInfo, info] = findInstruments(headers,template,info)

% extract variable names from templateStruct
variables = fields(template);

% loop through headers
for i = 1:length(headers)
    
    % find current header
    currentHeader = headers{i}(1,:);
    
    % loop through all variable names for template matches
    for j = 1:length(variables)
        
        % find current template in the loop
        currentTemplate = template.(variables{j});
        
        % check to make sure the template is not []
        if ~isempty(currentTemplate)
            
            templateLength = numel(currentTemplate);
            
            % find matches for current header and template
            currentMatches = strncmpi(currentTemplate(1:templateLength),currentHeader,templateLength);
            numMatches = sum(currentMatches);
            
            % find associated match heights
            temp = currentHeader(currentMatches);
            currentMatchHeights = nan(1,numMatches);
            for k = 1:numMatches
                localHeight = sscanf(temp{k}(templateLength+1:end),'%f');
                if ~isempty(localHeight)
                    currentMatchHeights(k) = sscanf(temp{k},[currentTemplate,'%f']);
                end
            end
            
            % if matches exist store them, along with height, in sensorInfo.  Column 1 is table number, Column 2 is
            % column number with in the table and column 3 is the height.
            if numMatches > 0
                if ~exist('sensorInfo','var')
                    sensorInfo = struct;
                end
                if isfield(sensorInfo,sprintf('%s',variables{j}))
                    sensorInfo.(variables{j})(end+1:end+numMatches,1) = i;
                    sensorInfo.(variables{j})(end-numMatches+1:end,2) = find(currentMatches == 1);
                    sensorInfo.(variables{j})(end-numMatches+1:end,3) = currentMatchHeights;
                else
                    sensorInfo.(variables{j})(1:numMatches,1) = i;
                    sensorInfo.(variables{j})(1:numMatches,2) = find(currentMatches == 1);
                    sensorInfo.(variables{j})(1:numMatches,3) = currentMatchHeights;
                end
            end
        end
    end
end

if ~exist('sensorInfo','var')
    sensorInfo = struct();
end
% output number of found instruments
for i = 1:length(variables)
    variable = variables{i};
    if isfield(sensorInfo,variable)
        numSensors = size(sensorInfo.(variables{i}),1);
    else
        numSensors = 0;
    end
    display(sprintf('variable name: %s.  Sensors found: %g',variable,numSensors));
end

flag1 = input('\nIs the number of discovered sensors correct(1 for yes, 0 for no): ');
if ~flag1
    error('Check templates!')
end
clc
% output sonic details
if isfield(sensorInfo,'u')
    numSonics = size(sensorInfo.u,1);
    for i = 1:numSonics
        height = sensorInfo.u(i,3);
        orientation = info.sonicOrientation(i);
        if info.sonicManufact(i)
            manufact = 'Campbell';
        else
            manufact = 'RMYoung ';
        end
        sensorInfo.u(i,4) = orientation;
        sensorInfo.u(i,5) = info.sonicManufact(i);
        display(sprintf('Sonic%g: Manufacturer=%s -  height=%gm - Orientation=%gdeg',i,manufact,height,orientation))
    end
    flag2 = input('\nIs the Sonic Information Correct (1 for yes, 0 for no): ');
    if ~ flag2
        error('Program stopped by user.  Check header format and height information')
    end
    % store sonic height in info struct
    info.sonicHeight = sensorInfo.u(:,3)';
else
    flag3 = input('\nNo sonics found.  Is this correct? (1 for yes, 0 for no): ');
    if ~ flag3
        error('Program stopped by user.  Check header format and height information')
    end
end

info = orderfields(info);
clc
end

