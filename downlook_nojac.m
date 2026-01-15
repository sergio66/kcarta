function [rads] = matlab_kcarta_downlook_nojac(head,prof,aux_struct,ropt0);
%
% function [rads]=matlab_kcarta_downlook_nojac(head,prof,aux_struct,ropt0);
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
%                    also whether to do SARTA NLTE or FAST kComp NLTE
% output
%   rads.freqAllChunks                  1x20000              freq
%   rads.radAllChunks                   20000x1              radiances
%   rads.iaa_kcomprstats_AllChunks          2x73             Sing Vectors stats

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

% edit this list to only keep gases you DO want!
% gasids = 1
% disp('>>>>>>>>>>>>>>>> only using gases ...')
% gasids = 103

fchunk = fA : df*10000 : fB; nchunk = length(fchunk);
absc = []; zang = [];

freqAllChunks = zeros(1,10000*nchunk);
radAllChunks  = zeros(10000*nchunk,1);
if iDebug > 0
  abscAllChunks   = zeros(10000*nchunk,nlays);
  abscAllChunksG2 = zeros(10000*nchunk,nlays);
  end

iaa_kcomprstats_AllChunks = [];

profileG = op_rtp_to_lbl2(1, refpro.glist, head, prof, refpro);

%% determine temp interp weights for each compressed chunk
[itlo,ithi,twlo,twhi,pi1Out] = temp_interp_weights(head,prof,refp);

%% determine temp interp weights for each continuum chunk
% freqX = (1:10000);  freqX = f0 + (freqX-1)*df;
% copt.cvers = CKD;
% copt.cdir  = cdir;
% copt.cswt  = cswt;
% copt.cfwt  = cfwt;
% [ci1,ci2,ctw1,ctw2] = continuum_temp_interp_weights(profileG, freqX, copt);
freqX = (1:10000);  freqX = f0 + (freqX-1)*df;
UseCKDTemp  = '5';
iUseCKDTemp = -1;
copt.cvers = CKD;
if str2num(CKD) < 0
  iUseCKDemp = +1;
  copt.cvers = UseCKDTemp;
end
copt.cdir  = cdir;
copt.cswt  = cswt;
copt.cfwt  = cfwt;
[ci1,ci2,ctw1,ctw2,cjtwlo,cjtwhi] = ...
   continuum_temp_interp_weights_jac(profileG, freqX, copt);
if str2num(CKD) < 0
  copt.cvers = CKD;
end

%tic
%profile on -history

for cc = 1 : length(fchunk)
  iaCountNumVec = [];

  ff = fchunk(cc);
  fr0 = ff + (0:9999)*df;
  fprintf(1,'doing chunk %4i \n',ff(1));

  absc = zeros(10000,nlays);

  chunkindex = (1:10000) + (cc-1)*10000;

  iDidNLTE = -1;
  for jj = 1 : length(gasids)
    gid = gasids(jj);
    %fprintf(1,'chunk %4i, doing gasID = %3i \n',ff,gid);

    if iDebug > 0 & gid == 2
      abscG2 = zeros(10000,nlays);
    end

    if gid == 1
      [absc,freq,iNumVec] = kcmix2(itlo, ithi, twlo, twhi, pi1Out, gid, ...
                                 profileG,ff,ropt0,refp,fr0,absc, prefix);

      iaCountNumVec(jj) = iNumVec;

      if str2num(copt.cvers) > 0
        cont = contcalc2(profileG,freq,copt,-1,ci1,ci2,ctw1,ctw2);
        absc = absc + cont; % disp(' >>>>>>>>>>>>> Adding cont!!! <<<<<<<<<<<')
        % absc = cont; disp(' >>>>>>>>>>>>> ONLY cont!!! <<<<<<<<<<<')
        % absc = absc; disp(' >>>>>>>>>>>>> NO cont!!! <<<<<<<<<<<')
      else
	disp('CKD = -1 so no continuum')
	cont = 0 * absc;	
      end	

    elseif gid == 2 & prof.solzen < 90 & ropt0.iNLTE == -2 & ...
       (ff(1) >= 2205 & ff(end) < 2405)

      fprintf(1,'CO2 FAST Compressed NLTE for chunk %4i solzen %8.6f iNLTE %2i\n',ff(1),prof.solzen,ropt0.iNLTE)

      iDidNLTE = +1;
      [abscOUT,freq,iNumVec] = kcmix2_NLTE(itlo, ithi, twlo, twhi, pi1Out, ...
                             gid,profileG,ff,ropt0,refp,fr0,absc, prefix,1);
      iaCountNumVec(jj) = iNumVec;
      if iDebug > 0
        abscG2 = abscOUT;
      end
      absc = absc + abscOUT;
%      semilogx(abscOUT',1:nlays)
%      ret

      [rplanckmod,freq,iNumVec] = kcmix2_NLTE(itlo, ithi, twlo, twhi, pi1Out,...
                           gid,profileG,ff,ropt0,refp,fr0,absc, prefix,2);
      %pcolor(freq,1:nlays,rplanckmod'); shading flat; colorbar;
%      semilogx(rplanckmod',1:nlays)
%      ret

    else
      %% CO2 LTE as well as all other gases 3-81
      [absc,freq,iNumVec] = kcmix2(itlo, ithi, twlo, twhi, pi1Out, gid, ...
                                 profileG,ff,ropt0,refp,fr0,absc, prefix);

      iaCountNumVec(jj) = iNumVec;
      
      end  
        
    end    %%% loop on gases

  absc(:,nlays) = absc(:,nlays)*rFracBot;
  if iDidNLTE ~= +1
    rplanckmod = ones(size(absc));
  end

  if iDebug > 0
    abscAllChunks(chunkindex,:) = absc;  
    abscAllChunksG2(chunkindex,:) = abscG2;  
    planckAllChunks(chunkindex,:) = rplanckmod;  
    end

  profX = prof;
  if iDebug > 1000 & iDebug < 2000
    profX.ptemp(ind_perturb) = profX.ptemp(ind_perturb) + 1;
    [prof.ptemp(ind_perturb) profX.ptemp(ind_perturb)]
    end

  [rad25,therm25,zang] = rtchunk_Tsurf(profX, absc, freq , rplanckmod, ropt0);

  freqAllChunks(chunkindex) = freq;
  radAllChunks(chunkindex) = rad25;

  iaa_kcomprstats_AllChunks = [iaa_kcomprstats_AllChunks; iaCountNumVec];
  end

%profile off
%toc
%keyboard


%%%%%%%%%%%%%%%%%%%% uncompression guts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear jacQGAllChunks jacTGAllChunks 
clear k2AllChunks jacQGAllChunks jacTGAllChunks

rads.freqAllChunks             = freqAllChunks;
rads.radAllChunks              = radAllChunks;
rads.gaslist                   = gasids;    
if iDebug > 0
  rads.abscAllChunks             = abscAllChunks;
  rads.abscAllChunksG2           = abscAllChunksG2;
  rads.planckAllChunks           = planckAllChunks;
  end
rads.iaa_kcomprstats_AllChunks = iaa_kcomprstats_AllChunks;

