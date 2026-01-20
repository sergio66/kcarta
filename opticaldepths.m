function ods = matlab_kcarta_opticaldepths(head,prof,aux_struct,ropt0,iGasDoOD,iBreakoutCont);
%
% fcn [absc]=matlab_kcarta_opticaldepths(head,prof,aux_struct,ropt0,iGasDoOD,iBreakoutCont);
%
% input
%   head,prof = head,prof of profile to be processed
%
%   aux_struct : set in set_aux
%     aux_struct.atm.nlays            number of layers
%     aux_struct.atm.rFracBot         fraction of bottom layer
%     aux_struct.atm.fA,fb,prefix     start and stop kCARTA chunks, prefix
%     aux_struct.cont.CKD,cswt,cfwt   continuum version, self/forn wgts
%     aux_struct.refp                 file with structure ref profiles, 
%                                     list of gasids
%   ropt0        = structure set in set_dirs, 
%                    with kpath,soldir,cdir,nltedir,co2chiDir
%                    also whether or not to do Thermal BackGnd and solar
%  iGasDoOD      = list of gases to do (or not to do)
%  iBreakoutCont = should we separately provide self,forn WV continuum components?
% 
% output
%   "y" is number of kCARTA chunks, "N" is number of layers
%   ods.freqAllChunks                  1xy0000              freq
%   ods.iaa_kcomprstats_AllChunks          2x73             Sing Vectors stats
%   ods.abscTotal                      y0000xN              optical depths, ALL gases
%   if giasID includes 1 & iBreakoutCont == 1
%     ods.abscS                        y0000xN              optical depths, self continuum
%     ods.abscF                        y0000xN              optical depths, forn continuum

%%%%%%%%%%%%%%%%%%%% uncompression guts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% plus clear sky rad transfer %%%%%%%%%%%%%%%%%%%%%%%%%%%%
co2ChiFilePath = ropt0.co2ChiFilePath;
cdir           = ropt0.cdir;

prefix   = aux_struct.atm.prefix;
fA       = aux_struct.atm.fA;
fB       = aux_struct.atm.fB;
df       = aux_struct.atm.df;
f0       = aux_struct.atm.f0;
nlays    = aux_struct.atm.nlays;
rFracBot = aux_struct.atm.rFracBot;
CKD      = aux_struct.cont.CKD;
cswt     = aux_struct.cont.cswt;
cfwt     = aux_struct.cont.cfwt;
refp     = aux_struct.refp;

global iDebug;
ind_perturb = -1;
if iDebug > 1000 & iDebug < 2000
  ind_perturb = iDebug - 1000;    %% perturb LAYTMP this layer if iDebug > 1000
elseif iDebug > 2000 & iDebug < 3000
  ind_perturb = iDebug - 2000;    %% perturb GASAMT this layer if iDebug > 2000
end
iGasDebug = 2;

%gasids = [[1 : 29],[51 : 63]];
load(refp);
gasids = refpro.glist;

if iGasDoOD == 9999
  %% no need to do anything
  disp('using ALL gases found in every chunk')
elseif iGasDoOD(1) > 0 & iGasDoOD(1) < 105
  %% use only these gases
  disp('using only gases found in list : ');
  abs(iGasDoOD)
  gasids = abs(iGasDoOD);
elseif iGasDoOD(1) < 0 
  %% use all EXCPET these gases
  disp('using all gases EXCEPT those found in list : ');
  abs(iGasDoOD);
  gasidsX = abs(iGasDoOD);
  gasids  = setdiff(gasids,gasidsX)
end

% edit this list to only keep gases you DO want!
% gasids = 1
% disp('>>>>>>>>>>>>>>>> only using gases ...')
% gasids = 103

fchunk = fA : df*10000 : fB; nchunk = length(fchunk);
absc = []; zang = [];

freqAllChunks = zeros(1,10000*nchunk);
abscAllChunks = zeros(10000*nchunk,nlays);
if length(intersect(1,gasids) == 1) & (iBreakoutCont == +1)
  selfAllChunks = zeros(10000*nchunk,nlays);
  fornAllChunks = zeros(10000*nchunk,nlays);
end

iaa_kcomprstats_AllChunks = [];

profileG = op_rtp_to_lbl2(1, refpro.glist, head, prof, refpro);

%% determine temp interp weights for each compressed chunk
[itlo,ithi,twlo,twhi,pi1Out] = temp_interp_weights(head,prof,refp);

%% determine temp interp weights for each continuum chunk
freqX = (1:10000);  freqX = f0 + (freqX-1)*df;
CKD0 = CKD;
if str2num(CKD) < 0
  CKD = '1';
end  
copt.cvers = CKD;
copt.cdir  = cdir;
copt.cswt  = cswt;
copt.cfwt  = cfwt;
[ci1,ci2,ctw1,ctw2] = continuum_temp_interp_weights(profileG, freqX, copt);
CKD = CKD0;
copt.cvers = CKD;

%tic
%profile on -history

for cc = 1 : length(fchunk)
  iaCountNumVec = [];

  ff = fchunk(cc);
  fr0 = ff + (0:9999)*df;
  fprintf(1,'doing chunk %4i \n',ff(1));

  absc = zeros(10000,nlays);

  chunkindex = (1:10000) + (cc-1)*10000;

  for jj = 1 : length(gasids)
    gid = gasids(jj);
    %fprintf(1,'   chunk %4i, doing gasID = %3i \n',ff,gid);

    [absc,freq,iNumVec] = kcmix2(itlo, ithi, twlo, twhi, pi1Out, gid, ...
                                 profileG,ff,ropt0,refp,fr0,absc,prefix);

    iaCountNumVec(jj) = iNumVec;

    if gid == 1 & str2num(copt.cvers) > 0
      if iBreakoutCont < 0
        cont = contcalc2(profileG,freq,copt,ci1,ci2,ctw1,ctw2);
      elseif iBreakoutCont == +1
        [cont,contS,contF]=contcalc2_S_F(profileG,freq,copt,ci1,ci2,ctw1,ctw2);
      end
      absc = absc + cont; disp(' >>>>>>>>>>>>> Adding cont!!! <<<<<<<<<<<')
      % absc = cont; disp(' >>>>>>>>>>>>> ONLY cont!!! <<<<<<<<<<<')
      % absc = absc; disp(' >>>>>>>>>>>>> NO cont!!! <<<<<<<<<<<')
    end

  end    %%% loop on gases

  % absc(:,nlays) = absc(:,nlays)*rFracBot;
  freqAllChunks(chunkindex) = freq;
  abscAllChunks(chunkindex,:) = absc;  
  if iBreakoutCont == 1 & length(intersect(1,gasids)== 1)
    selfAllChunks(chunkindex,:) = contS;  
    fornAllChunks(chunkindex,:) = contF;  
  end

  iaa_kcomprstats_AllChunks = [iaa_kcomprstats_AllChunks; iaCountNumVec];
end

%profile off
%toc
%keyboard


%%%%%%%%%%%%%%%%%%%% uncompression guts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear jacQGAllChunks jacTGAllChunks iGasDoOD
clear k2AllChunks jacQGAllChunks jacTGAllChunks

ods.gaslist                   = gasids;
ods.freqAllChunks             = freqAllChunks;
ods.abscTotalAllChunks        = fliplr(abscAllChunks);
if iBreakoutCont == 1 & length(intersect(1,gasids)== 1)
  ods.selfAllChunks           = fliplr(selfAllChunks);
  ods.fornAllChunks           = fliplr(fornAllChunks);
end
ods.iaa_kcomprstats_AllChunks = iaa_kcomprstats_AllChunks;

disp(' ')
disp('did fliplr of abscAllChunks, selfAllChunks,fornAllChunks')
