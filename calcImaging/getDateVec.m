function [Y M D H MN S] = getDateVec(datestr)
    %This is a function for extracting recording time info from the
    %Finchscope filenames
    Y = (datestr(1:4));
    M = (datestr(6:7));
    D = (datestr(9:10));
    H = (datestr(12:13));
    MN = (datestr(15:16));
    S = (datestr(18:19));
end