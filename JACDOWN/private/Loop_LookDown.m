% this subroutine does the general looping for the Jacobians, so all that
% has to be called is raaResults with the appropriate raaDT or raaaDQ
function raaGeneral = Loop_LookDown(rSatAngle,raLayAngles,raSunRefl,...
             raUseEmissivity,raSurface,raSun,raSunAngles,raThermal,...
             raaOneMinusTau,raaLay2Sp,raaRad,prof)

[aa,iNumLayer] = size(raaOneMinusTau);

rCos  = cos(prof.satzen*pi/180.0);
rWsun = cos(prof.solzen*pi/180.0);

% do the bottommost layer first
iLyr = 1;
% first do the surface term
iJ1 = 1;
raTemp = JacobTerm(iJ1,iLyr,raaLay2Sp);
raaGeneral(:,iLyr) = raUseEmissivity.*raSurface.*raTemp;

% recall raTemp is raL2S from gnd to top
iLay = 1;
rCos  = cos(raLayAngles(iLay)*pi/180.0);
rWsun = cos(raSunAngles(iLay)*pi/180.0);

kSolar = -1;
if raSunAngles(iLay) < 90
  kSolar = 1;
  raaGeneral(:,iLyr) = raaGeneral(:,iLyr) + ...
		       raSunRefl.*raSun.*raTemp.*(1+rCos/rWsun);
end

% include the EASY part of thermal contribution
kThermal = -1;
if raThermal(1) > 0
  kThermal = +1;
  raaGeneral(:,iLyr) = raaGeneral(:,iLyr) + ...
		       (1.0-raUseEmissivity)/pi.*raThermal.*raTemp;
end

% set raTemp1 to all zeros (since this is the bottom layer, there is no 
% cumulative contribution
raTemp1 = zeros(size(raTemp));

% loop over the remaining layers
for iLyr = 2 : iNumLayer
  iLay = iLyr;
  rCos = cos(raLayAngles(iLay)*pi/180.0);
  % first do the surface term
  iJ1 = 1;
  raTemp = JacobTerm(iJ1,iLyr,raaLay2Sp);
  raaGeneral(:,iLyr) = raUseEmissivity.*raSurface.*raTemp;
  
  % recall raTemp is raL2S from gnd to top
  if (kSolar > 0) 
    rWsun=cos(raSunAngles(iLay)*pi/180.0);
    raaGeneral(:,iLyr) = raaGeneral(:,iLyr) + ...
      raSunRefl.*raSun.*raTemp*(1+rCos/rWsun);
    end

  % include the EASY part of thermal contribution
  if (kThermal > 0) 
    raaGeneral(:,iLyr) = raaGeneral(:,iLyr) + ...
      (1.0-raUseEmissivity)/pi.*raThermal.*raTemp;
    end

  % now loop over the layers that contribute (i.e. < iLyr) ....
  iJ = iLyr - 1;
  iJ1 = iJ + 1;
  raTemp = JacobTerm(iJ1,iLyr,raaLay2Sp);
  raTemp1 = raTemp1 + raaOneMinusTau(:,iJ).*raaRad(:,iJ).*raTemp;
  raaGeneral(:,iLyr) = raaGeneral(:,iLyr) + raTemp1;
  
end
