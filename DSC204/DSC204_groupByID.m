function dscgroups = DSC204_groupByID(dscdata)
% dscgroups = DSC204_groupByID(dscdata)
%
% Groups the DSC204 data in the given structure array by their ID in individual cells.
%
% INPUT:    dscdata --> DSC data as read by DSC204_readFile(s)
%   
% OUTPUT: dscgroups --> cell array, with each cell containing the dsc data 
%                       belonging to identical ID
%
%
% Author: Andreas Sommer, Apr2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%


% get the IDs
IDs = arrayfun(@(x) x.desc.IDENTITY, dscdata, 'UniformOutput', false);

% get unique ids, sort them
IDs = sort(unique(IDs));

% now assemble the groups
dscgroups = cell(length(IDs), 1);
for k = 1:length(dscdata)
   idx = find(strcmp(dscdata(k).desc.IDENTITY, IDs));
   dscgroups{idx}(end+1) = dscdata(k); %#ok<FNDSB> It's fast enough here
end



end