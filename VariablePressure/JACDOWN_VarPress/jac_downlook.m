function [qjac,tjac,wgt,sjac,ejac] = ...
  jac_downlook(raFreq,zang,efine,raSol0,raThermal,absc,raaRad,...
               jacQG,jacTG,prof,iDoJac);

% computes jacs and weighting functions
% input
%   efreq     = 1e4 x 1          input wavenumbers
%   efine     = 1e4 x 1          input emissivity
%   raSol     = 1e4 x 1          input solar radiance at gnd
%   raThermal = 1e4 x 1          input thermal backgnd radiance at gnd
%   zang   = nlays x 1           input layer angles
%   absc   = 1e4 x nlays         gas ODs
%   jacQG  = 1e4 x nlays         dK/dq for needed gas
%   jacTG  = 1e4 x nlays         dK/dT from all gases
%   prof                         structure with info such as : temp etc
%   
% output 
%    qjac = 1e4 x nlays      dr/dq
%    tjac = 1e4 x nlays      dr/dT
%    wgt  = 1e4 x nlays      Wgt
%    sjac = 1e4 x 1          dr/d(SurfTemp)
%    ejac = 1e4 x 1          dr/d(SurfEmis)

global iDebug

[aa,iNumLayer] = size(absc);
if iNumLayer ~= prof.nlevs-1
  error('inconsistent number of layers')
  end

if prof.plevs(1) < prof.plevs(10)
  %% pressures increasing with index ==> height decreasing with index
  %% ==> index 1 = TOA
  %% f77 kCARTA assumed index 1 = GND, so flip
  raVtemp   = flipud(prof.ptemp(1:prof.nlevs-1));
  absc      = fliplr(absc);
  jx = jacQG;
  for ii = 1 : length(iDoJac)
    jox = squeeze(jx(ii,:,:));
    jox = fliplr(jox);
    jacQG(ii,:,:) = jox;
    end
  jacTG     = fliplr(jacTG);
  rTSurface = prof.stemp;
  %zang      = flipud(zang);
else
  %% no need to flip
  raVtemp   = (prof.ptemp(1:prof.nlevs-1));
  absc      = (absc);
  jacQG     = (jacQG);
  jacTG     = (jacTG);
  rTSurface = prof.stemp;
  zang      = (zang);
  end

%% initialise the layer-to-space matrix
raaLay2Sp = AtmosLayer2Space(absc,zang);
if iDebug == 1
  for ii = iNumLayer : -1 : 1
    fprintf(1,'%3i %8.4f %8.6f  %8.6f \n',ii,zang(ii),absc(1,ii),raaLay2Sp(1,ii))
    end
  end

%disp('initializing Jac radiances/d/dT(radiances) ...')
[raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd] = ...
  DoPlanck_LookDown(prof,raFreq,zang,absc,raVtemp,raaRad);

%disp('initializing Jacobian loops ...')
raSunAngles = vaconv(prof.solzen,prof.zobs,prof.palts);
raSurface = ttorad(raFreq,prof.stemp);
raaGeneral = Loop_LookDown(prof.satzen,zang,(1-efine)/pi,efine,...
             ttorad(raFreq,prof.stemp),raSol0,raSunAngles,raThermal,...
             raaOneMinusTau,raaLay2Sp,raaRad,prof);

kThermal = +1;
kThermalJacob = +1;
raaGeneralTh = zeros(size(raaGeneral));
raaOneMinusTauTh = zeros(size(raaOneMinusTau));
if ((kThermal >=  0) & (kThermalJacob >  0))
  %disp('initializing Jacobian thermal loops ...')
  [raaOneMinusTauTh,raaGeneralTh] = Loop_thermaldown(prof.satzen,zang,...
                                          absc,raFreq,raaRad,raaLay2Gnd);
  end

for ii = 1 : length(iDoJac)
  if iDebug == 1
    fprintf(1,'\n gas  d/dq : jacID%2i \n',ii)
  else
    fprintf(1,'\n gas  d/dq : ')
    end
  for iLay=1 : iNumLayer
    rWeight = 1;
    fprintf(1,'.')
    qjac(ii,:,iLay) = JacobGasAmtFM1(iLay,raFreq,raaRad,raaRadDT,efine,...
                        raaOneMinusTau,raaTau,squeeze(jacQG(ii,:,:)),...
                        raaLay2Sp,raThermal,raaLay2Gnd,...
                        prof.satzen,zang,...
                        raaGeneral,raaGeneralTh,raaOneMinusTauTh);
     end
   end

if iDebug == 1
  fprintf(1,'\n temp d/dT : \n')
else
  fprintf(1,'\n temp d/dT : ')
  end
for iLay=1 : iNumLayer
  fprintf(1,'.')
  rWeight = 1.0;
  tjac(:,iLay) = JacobTempFM1(iLay,raFreq,raaRad,raaRadDT,efine,...
                   raaOneMinusTau,raaTau,jacTG,...
                   raaLay2Sp,raThermal,raaLay2Gnd,...
                   prof.satzen,zang,...
                   raaGeneral,raaGeneralTh,raaOneMinusTauTh);
   end

fprintf(1,'\n wgt       : ')
for iLay=1 : iNumLayer
  fprintf(1,'.')
  rWeight = 1.0;
  wgt(:,iLay) = wgtfcndown(iLay,prof.satzen,zang,raaLay2Sp,absc);
  end

sjac = surface_temp_jacobian(raFreq,prof.stemp,efine,raaLay2Sp);
ejac = surface_emis_jacobian(raFreq,prof.stemp,efine,raThermal,raaLay2Sp);

fprintf(1,' \n');