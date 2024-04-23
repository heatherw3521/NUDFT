% Set up MATLAB environment and paths.

% Barnett 11/4/22.

% for INUDFT, 9/14/23.
addpath .

% these apparently needed
addpath modifed_HMtoolbox_files2019
addpath extras

addpath tests
addpath itersolv  % where AHB working

% dependencies:
if ~exist('finufft1d2')
  warning('please set up your FINUFFT MATLAB path!');
  %addpath finufft/matlab     % try to add: put your paths to FINUFFT here!
end
