function header = importHeader(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   BEACH = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   BEACH = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   beach = importfile('beach.dat', 2, 2);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2013/11/26 12:27:00

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 5;
    endRow = 5;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
temp = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(temp)
        temp{col} = [temp{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);
header = cell(1,length(temp));

for i = 1:length(temp)
    if ~isempty(temp{i})
        header(i) = strrep(temp{i},'"','');
    end
end