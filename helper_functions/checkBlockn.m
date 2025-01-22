function subdirectories = checkBlockn(directory_path)

% Get the list of all contents (files and directories) in the directory
contents = dir(directory_path);

% Initialize a cell array to store the names of subdirectories
subdirectories = {};

% Loop through each item in the contents
for i = 1:numel(contents)
    % Check if the item is a directory and its name starts with "Block"
    if contents(i).isdir && startsWith(contents(i).name, 'Block')
        % If so, add its name to the subdirectories cell array
        subdirectories{end+1} = contents(i).name;
    end
end

end