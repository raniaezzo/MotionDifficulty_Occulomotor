function SamplingRate = findSamplingRate(filename)

    % find lines that contain "RATE" in msg file
    A = regexp(fileread(filename),'\n','split');
    whichline = find(contains(A,'RATE'));
    RATELINE = split(char(A(whichline(1))));
    SamplingRate = str2double(RATELINE{5});

end