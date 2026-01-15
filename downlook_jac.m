function [rads,jacs] = ...
    matlab_kcarta_downlook_jac(head,prof,aux_struct,ropt0,iDoJac,iJacobOutput);
% function 
%   [rads,jacs] = ...
%  matlab_kcarta_downlook_jac(head,prof,aux_struct,ropt0,iDoJac,iJacobOutput)
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
%   iDoJac       = list of gasIDS for jacs
%     iDoJac = -1;                      %% no jacs
%     iDoJac = [iGid1 iGid2 ... iGidN]; %% do jacobians; the iGid tells which
%                                       %% gases to do amt jacs for
%                                       %% (temp jacs always done)
%     warning : iGid1 = WV includes lines and continuum
%     example : iDoJac = [1 3];
%   iJacobOutput = what kind of jacobian output
%     iJacobOutput = -1;        %% dr/dT, dr/dq
%     iJacobOutput =  0;        %% dr/dT, dr/dq*q
%     iJacobOutput = +1;        %% dBT/dT, dBT/dq*q
%
% output
%   rads.freqAllChunks                  1x20000              freq
%   rads.radAllChunks                   20000x1              radiances
%   rads.iaa_kcomprstats_AllChunks          2x73             Sing Vectors stats
%     jacs.ejacAllChunks                  20000x1              emissivity jac
%     jacs.qjacAllChunks                  2x20000x96           gas jacs
%     jacs.sjacAllChunks                  20000x1              stemp jac
%     jacs.tjacAllChunks                  20000x96             tempr jacs
%     jacs.wgtAllChunks                   20000x96             wgt fcns
%

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
  abscAllChunks = zeros(10000*nchunk,nlays);
  end

%% q = gas amt, t = temp, s = surfacetemp, e = surface emiss
qjacAllChunks = zeros(length(iDoJac),nchunk*10000,nlays);
tjacAllChunks = zeros(nchunk*10000,nlays);
wgtAllChunks  = zeros(nchunk*10000,nlays);
sjacAllChunks = zeros(10000*nchunk,1);
ejacAllChunks = zeros(10000*nchunk,1);

iaa_kcomprstats_AllChunks = [];

profileG = op_rtp_to_lbl2(1, refpro.glist, head, prof, refpro);

%% determine temp interp weights for each chunk
[itlo,ithi,twlo,twhi,jtwlo,jtwhi,pi1Out] = ...
   temp_interp_weights_jac(head,prof,refp);

%% determine temp interp weights for each continuum chunk
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

for cc = 1 : nchunk
  iaCountNumVec = [];

  ff = fchunk(cc);
  fr0 = ff + (0:9999)*df;
  fprintf(1,'doing chunk %4i \n',ff(1));

  absc = zeros(10000,nlays);

  jacTchunk = zeros(10000,nlays);

  chunkindex = (1:10000) + (cc-1)*10000;

  for jj = 1 : length(gasids)
    gid = gasids(jj);
    %fprintf(1,'chunk %4i, doing gasID = %3i \n',ff,gid);

    if iDebug > 1000 & iDebug < 2000
      fprintf(1,'perturb T( %3i ) by 1 K \n',ind_perturb);
      profileG.mtemp(ind_perturb) = profileG.mtemp(ind_perturb) + 1;
      profileG.mtemp(ind_perturb)
    elseif iDebug > 2000 & gid == iGasDebug
      fprintf(1,'perturb Q %2i (%3i) by 10 percent \n',iGasDebug,ind_perturb);
      profileG.gamnt(ind_perturb)*0.1
      profileG.gamnt(ind_perturb) = profileG.gamnt(ind_perturb)*1.1;
      end

    iGasExist = -1;
    [absc,freq,iNumVec,gasprofX,jacTG,jacQG,iGasExist] = ...
               kcmix2jac(itlo,ithi,twlo,twhi,jtwlo,jtwhi,pi1Out,gid,iDoJac, ...
                              profileG,ff,ropt0,refp,fr0,absc,prefix);

    [aa,iNumLayer] = size(absc);
    if length(intersect(iDoJac,gasids(jj))) == 1
      iQG = find(gasids(jj) == iDoJac);
      gasprofQG(iQG,:,:) = gasprofX;
      end

    iaCountNumVec(jj) = iNumVec;

    if gid == 1
      if str2num(copt.cvers) > 0	
        [cont,contjacT,contjacQ] = ...
         contjaccalc2(profileG,freq,copt,iDoJac,ci1,ci2,ctw1,ctw2,...
                      cjtwlo,cjtwhi);
        jacTG = jacTG + contjacT;  
        if length(intersect(iDoJac,1)) == +1
          jacQG = jacQG + contjacQ; %% water lines + water cont jac
        end
        absc = absc + cont;
      else
	disp('CKD = -1 so no continuum')
	cont = 0 * absc;	
	contjacT = 0 * jacTG;
        if length(intersect(iDoJac,1)) == +1
          contjacQ = 0 * jacQG; %% water lines + water cont jac
        end          	
      end
    end
      
    if iGasExist == +1
      jacTchunk = jacTchunk + jacTG;
      end

    if length(intersect(iDoJac,gasids(jj))) == 1
      iQG = find(gasids(jj) == iDoJac);
      jacQchunk(iQG,:,:) = jacQG;
    end    

  end    %%% loop on gases

  absc(:,nlays) = absc(:,nlays)*rFracBot;
  if iDebug > 0
    abscAllChunks(chunkindex,:) = absc;  
  end

  [rad25,therm25,zang,efine,rSol25,raaRad] = ...
     rtchunk_Tsurf_jac(prof, absc, freq , ropt0);
  [qjac,tjac,wgt,sjac,ejac] = ...
       jac_downlook(freq,zang,efine,rSol25,therm25,absc,raaRad,...
                      jacQchunk,jacTchunk,prof,iDoJac);
  if iJacobOutput == 0
    [aa,iNumLayer] = size(absc);
    qjac = doQjacOutput(gasprofQG,qjac,iDoJac,0,iNumLayer,freq,rad25);
  elseif iJacobOutput == 1
    [aa,iNumLayer] = size(absc);
    qjac = doQjacOutput(gasprofQG,qjac,iDoJac,1,iNumLayer,freq,rad25);
    tjac = tjac .* (dbtdr(freq,rad25')'*ones(1,iNumLayer));
    sjac = sjac .* (dbtdr(freq,rad25')');
    ejac = ejac .* (dbtdr(freq,rad25')');
  end

  qjacAllChunks(:,chunkindex,:) = qjac;
  tjacAllChunks(chunkindex,:)   = tjac;
  wgtAllChunks(chunkindex,:)    = wgt;
  sjacAllChunks(chunkindex)     = sjac;
  ejacAllChunks(chunkindex)     = ejac;

  freqAllChunks(chunkindex) = freq;
  radAllChunks(chunkindex) = rad25;

  iaa_kcomprstats_AllChunks = [iaa_kcomprstats_AllChunks; iaCountNumVec];
end
%%%%%%%%%%%%%%%%%%%% uncompression guts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%profile off
%toc
%keyboard

rads.freqAllChunks             = freqAllChunks;
rads.radAllChunks              = radAllChunks;
rads.gaslist                   = gasids;    
if iDebug > 0
  rads.abscAllChunks             = abscAllChunks;
  end
rads.iaa_kcomprstats_AllChunks = iaa_kcomprstats_AllChunks;

jacs.ejacAllChunks = ejacAllChunks;
jacs.qjacAllChunks = qjacAllChunks;
jacs.sjacAllChunks = sjacAllChunks;
jacs.tjacAllChunks = tjacAllChunks;
jacs.wgtAllChunks  = wgtAllChunks;
