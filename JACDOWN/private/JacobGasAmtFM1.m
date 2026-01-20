% this subroutine does the Jacobians wrt gas amount
function raResults = JacobGasAmtFM1(...
         iLay,raFreq,raaRad,raaRadDT,raUseEmissivity,...
         raaOneMinusTau,raaTau,jacQG,raaLay2Sp,...
         raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,...
         raaGeneralTh,raaOneMinusTauTh);

global iDebug

[aa,iNumLayer] = size(raaRad);

iM1 = iLay;

% fix the sat angle weight factor
rSec = 1.0/cos(rSatAngle*pi/180.0);
rSec = 1.0/cos(raLayAngles(iM1)*pi/180.0);

% read the appropriate layer from general results
% this includes all the surface terms 
raResults = raaGeneral(:,iLay);

% set the constant factor we have to multiply results with
% this is a gas amt jacobian
raTemp = jacQG(:,iM1);

raResults = MinusOne(raTemp,raResults);

% add on the the derivatives wrt radiating layer
if (iLay <  iNumLayer) 
  % this is not the topmost layer
  iJ1 = iLay;
  raResults = raResults + raTemp.*raaRad(:,iJ1).*raaLay2Sp(:,iJ1);
elseif (iLay == iNumLayer)
  % do the topmost layer correctly
  iJ1 = iLay;
  raResults = raResults + raTemp.*raaTau(:,iJ1).*raaRad(:,iJ1);
  end

% now multiply results by the 1/cos(viewing angle) factor
if (abs(rSec-1.0000000) >= 1.0E-5) 
  raResults = raResults*rSec;
  end

iDoBckGnd = +1;
if iDoBckGnd > 0
  raResultsTh = JacobThermalAmtFM1(raFreq,raaRad,iLay,...
                 raUseEmissivity,raTemp,raaLay2Sp,...
                 raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh);
  %plot(1:10000,raResults,1:10000,raResultsTh,'g')
  %pause(0.1)
  raResults = raResultsTh + raResults;
  end

iJunk = 1;
iJunk = 8402;
iDebugX = iDebug;
iDebug = +1;
iDebug = -1;
if iDebug == 1
  %      print *,iMMM,iM1,raaGeneral(iJunk,iMMM),raaAllDT(iJunk,iM1),raTemp(iJunk),raResults(iJunk),
  %     $         raaRad(iJunk,iMMM),raaRadDT(iJunk,iMMM),raaLay2Gnd(iJunk,iMMM),raaOneMinusTau(iJunk,iMMM)
  data = [-1 iLay raaGeneral(iJunk,iLay) jacQG(iJunk,iLay) raTemp(iJunk) ...
                 raaRad(iJunk,iLay) raaRadDT(iJunk,iLay) raaLay2Sp(iJunk,iLay) raaOneMinusTau(iJunk,iLay) ...
                 raResults(iJunk)];
  fprintf(1,' %3i %3i %10.6e %10.6e %10.6e %10.6e %10.6f %10.6f %10.6f %10.6e \n',data);
  end
iDebug = iDebugX;
