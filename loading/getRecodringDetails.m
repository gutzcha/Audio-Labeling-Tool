function [details] = getRecodringDetails(file,path,str)
%[details] = getRecodringDetails(file,path,str)
%file - file name with extention
%path - file path without file name
%str - file name without extention

%   Split the name and set the details


str = lower(str);

%Get rat number
ratNum = str(strfind(str,'rat')+3);

%Get day number
dayNum = str(strfind(str,'day')+3);

%Get paradigm
if contains(str,'chamber')
    paradigm = 'Chamber';
elseif contains(str,'free')
    paradigm = 'Free';
else
    paradigm='';
end

details.ratNum = str2double(ratNum);
details.paradigm = paradigm;
details.dayNum = str2double(dayNum);
details.fileName = file;
details.filePath = path;
