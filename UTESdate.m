function [data, dataInfo, info] = date(data, dataInfo, info)
display(sprintf('\nStoring serial date number in first column of each table and deleting columns 2 - 4'))
for i = 1:length(data)
    display(sprintf('Finding serial dates for table %g',i))
    % check to make sure table loaded
    if ~isempty(data{i})
        
        % grab date columns
        dateVec = data{i}(:,1:4);
        
        % find year. If experiment crosses from leap year to non-leap year or visa versa, tables need to be split at 12-31.
        yr = dateVec(1,1);
        if yr == 2004 || yr == 2008 || yr == 2012 || yr == 2016
            doy = [0 31 60 91 121 152 182 213 244 274 305 335 366];
        else
            doy = [0 31 59 90 120 151 181 212 243 273 304 334 365];
        end
        month = zeros(length(dateVec),1);
        day = zeros(length(dateVec),1);
        year = dateVec(:,1);        % convert from DOY to MM DD
        
        % find unique days
        doys = unique(dateVec(:,2));
        
        % iterate through unique days
        for j = 1:length(doys)
            dayPointer = find(dateVec(:,2) == doys(j));
            localMonth = find(doy >= abs(doys(j)),1,'first')-1;
            localDay = abs(doys(j)) - doy(localMonth);
            
            % place day and month in correct rows
            month(dayPointer) = localMonth;
            day(dayPointer) = localDay;
        end
        
        hour = floor(dateVec(:,3)*10^-2);
        min = dateVec(:,3)-hour*10^2;
        sec = dateVec(:,4);
        
        % create serial date numbers
        date = datenum(year,month,day,hour,min,sec);
        
        % store in col 1 of table
        data{i}(:,1) = date;
        
        % delete columns 2-4
        data{i}(:,2:4) = [];
        
        % store begin and end date in dateInfo
        row = find(cellfun(@isempty,dataInfo(:,i))==1,1,'first');
        try
            if isempty(row)
                dataInfo{end+1,i} = strcat('beg date:',datestr(date(1)));
                dataInfo{end+1,i} = strcat('end date:',datestr(date(end)));
            else
                dataInfo{row,i} = strcat('beg date:',datestr(date(1)));
                dataInfo{row+1,i} = strcat('end date:',datestr(date(end)));
            end
            
            info.date = datestr(date(1),'yyyy_mm_dd');
            
        catch err
            warning(err.message)
        end
    end
end
end

