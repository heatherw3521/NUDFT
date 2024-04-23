function xy = mergestructs(x,y)
% MERGESTRUCTS  combine two structs with distinct fields
%
% xy = mergestructs(x,y) returns single structs combining fields of two
%  structs. Raises error if any field names are common to the inputs.
%
% With no arguments does self-test
%
% See: https://www.mathworks.com/matlabcentral/answers/96973-how-can-i-concatenate-or-merge-two-structures
if nargin==0, test_mergestructs; return; end
  
xy = cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

%%%%%%%
function test_mergestructs
mergestructs(struct('x',1,'y',2),struct('z',3))
try
  mergestructs(struct('x',1,'y',2),struct('x',3))   % duplicate field x
catch ME
  disp(['correctly errored: ',ME.message])
end
