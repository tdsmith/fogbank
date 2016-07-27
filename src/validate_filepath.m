% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

function [cur_path] = validate_filepath(cur_path)
% check the number of inputs
if nargin ~= 1, return, end
% check that the cur_path variable is a char string
if ~isa(cur_path, 'char')
    error('validate_filepath:argChk','invalid input type');
end

% get the file attributes
[status,message] = fileattrib(cur_path);
% if status is 0 then the file path was invalid
if status == 0
    error('validate_filepath:notFoundInPath', 'No such file or directory: \"%s\"',cur_path);
else
	% cur_path held a valid file path to either a directory or a file
    cur_path = message.Name;
    % determine if cur_path is a file or a folder
    if message.directory == 0
        % the path represents a file
        % do nothing
    else
        % the path represents a directory
        if cur_path(end) ~= filesep
            cur_path = [cur_path filesep];
        end
    end
end

end