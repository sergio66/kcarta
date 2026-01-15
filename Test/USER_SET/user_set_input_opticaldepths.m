% this driver file is used to set up the input parameters for a 
% optical depthcomputation
% Assumption : user_set_dirs.m has been set correctly!!!!
%
% S. Machado, UMBC, March 2011        sergio@umbc.edu
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS (with examples)
%
% iAirs drives klayers, and eventually the convolver
%   iAirs = 0,-1,+1,2,3  : raw, AIRS, IASI, CRiS
%   example : iAirs = +1;     %%% doing things for AIRS
%
% iHITRAN sets the kCompressed directory, based on HITRAN
%   iHITRAN = 2000,2004,2008 are your options
%
% iMatlab_vs_f77 tells the code whether the database is Matlab or ieee-le 
%   iMatlab_vs_f77 = +1   use Matlab version
%   iMatlab_vs_f77 = -1   use f77 version
%
% fA and fB are start and stop wavenumbers
%   example : fA = 605;  fB = 830; 
%
% iGasDoOD controls the optical depths
%    +9999                      do all
%    [gid1   gid2 ...  gidN]    only do gasIDs gid1,gid2, ... gidN
%    [-gid1 -gid2 ... -gidN]    do all except gasIDs gid1,gid2, ... gidN
% warning gasID = 1 (WV) automatically includes continuum
% warning gasID = 2 (CO2) never includes NLTE OD contributions
%
% CKD is the CKD version : choose 1,2,3,4,5
%   example : CKD = '1';
%
% iBreakoutCont tells the code whether or not to output SLEF/FORN continuum as
% as well as TOTAL OD       -1 : no  (only output TOTAL)
%                           +1 : yes (output TOTAL, SELF, FORN)
%   example iBreakOutCont = +1
%
% these next two define where the input rtp file is in
%   dir is the directory where the file is in
%   fin is the actual file
%   example : dirin = '/asl/s1/sergio'; 
%             fin = 'pin_feb2002_sea_airsnadir_op.sun.rtp';
%
% iProfRun is which of the rtp profiles to run
%   example iProfRun = 1;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS (with examples)
%
% structures "odOut"  is your output variables
%            "ropt0"   and "stuff"   give some useful input info
%
% example if fA,fB = 605,805; each 10000 point kCARTA chunk is 25 cm-1 wide, 
% so there are 9 chunks or 90000 wavenumber points in this run
%
% odsOut : always produced
%    freqAllChunks                  1x90000          freq      cm-1
%    abscAllChunks                  90000xN          total OD for N layers
%    gasids                         1xG              list of gasIDS used
% iaa_kcomprstats_AllChunks          2x73            Singular Vectors stats
%
% ropt0 : reproduces important input parameters (from set_dirs) eg
%             kpath: '/asl/s1/sergio/KCMIX_DATABASE/H2004_matlab'
%            soldir: '/home/sergio/MATLABCODE/KCMIX2/KCMIXCODE/solarV2'
%              cdir: '/asl/data/kcarta/KCARTADATA/General/CKDieee_le'
%           nltedir: [1x60 char]
%    co2ChiFilePath: '/asl/data/kcarta/KCARTADATA/General/ChiFile/'
%            rsolar: 1
%            rtherm: 2
%
% stuff : not too important, duplicates stuff set here (set_input_downlook.m)
%            freqs: [605 805]
%    input_rtpfile: [1x66 char]
%       layersprof: [1x1 struct]
%         iGasDoOD: [1 2]
%         iProfRun: 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% user defined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% user defined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% user defined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iAirs drives klayers, and ultimately the convolver
iAirs = +1;

% iHITRAN sets the kCompressed directory, based on HITRAN
iHITRAN = 2008;

% iMatlab_vs_f77 : use Matlab (+1) or ieee-le (-1) kcomp database
%iMatlab_vs_f77 = +1;   % use Matlab version
iMatlab_vs_f77 = -1;   % use f77 version

% fA and fB are start and stop wavenumbers
%fA = 500   fB =  605; 
fA = 605;  fB = 2830; 
fA = 1780;  fB = 1805; 

% iGasDoOD gives the gasIDs to be included
%iGasDoOD = [1 103];        %% only do WV (includes continuum) and CO2
iGasDoOD = 9999;            %% do all
iGasDoOD = [9999];            %% do all

% CKD is the CKD version : choose 1,2,3,4,5
CKD = '6';

% iBreakoutCont tells whether to output SELF/FORN continuum as well
iBreakoutCont = +1;

% these next two define the input dir and rtp file
dirin = '/home/sergio/MATLABCODE/KCMIX2/PACKAGE_UPnDOWNLOOK_2011/RTPFILES/';
  fin   = 'desert_op.rtp'; 
  iProfRun = 1;
dirin = '/asl/s1/strow/rtprod_cris/2013/08/28/';
  fin   = 'test.rtp';
  iProfRun = 262;

