% this subroutine adds on the AIRS NLTE predictors. Note the AIRS channels
% are spaced approx 2 cm-1 apart in this region, so we need to interpolate
% to the 0.0025 cm-1 of KCARTA
function radNLTE = nlte(raFreq,satzen,raSatAngles,raSunAngles,raVT,iNumLayer,nltedir)

if raFreq(1) > 2405 | raFreq(length(raFreq)) < 2205
  %% no NLTE
  return
  end

suncos = raSunAngles(1);           %% at surface
scos1  = raSunAngles(iNumLayer);   %% at TOA
vsec1  = raSatAngles(iNumLayer);   %% at TOA

suncos = cos(suncos*pi/180.0);
scos1  = cos(scos1*pi/180.0);
vsec1  = 1/cos(vsec1*pi/180.0);

% this subroutine mimics the SARTA NLTE model
% see /asl/packages/sartaV106/Src/calnte.f

tHigh = sum(raVT(iNumLayer-4:iNumLayer))/5;

      pred(1) = 1.0;
      pred(2) = scos1;
      pred(3) = scos1*scos1;
      pred(4) = scos1*vsec1;
      pred(5) = scos1*tHigh;
      pred(6) = suncos;

%     !/asl/packages/sartaV106/Src/incFTC_apr05_nte.f
%     FNCOFN = '/asl/data/sarta_database/Data_jan04untun/Coef/setnte_oct05.dat'
%     !matlab cd : /home/sergio/KCARTADATA/NLTE/SARTA_COEFS/
%     !fname='/asl/data/sarta_database/Data_jan04untun/Coef/setnte_oct05.dat';
%     ![idchan, freq, coef] = rd_nte_be(fname);
%     !/wrt_nte_le.m
fncoefn = nltedir;
ee = exist(fncoefn);
if ee == 0
  fprintf(1,'you have set nlte file = %s \n',fncoefn)
  error('file DNE');
  end
[ichan, frqchn, coefn] = rd_nte_le(fncoefn);
nchnte = length(ichan);

% !/asl/packages/sartaV106/Src/calnte.f
raDrad = coefn * pred';

radNLTE = interp1(frqchn,raDrad,raFreq,'linear');
oo = find(raFreq > max(frqchn) | raFreq < min(frqchn));
if length(oo) > 0
  radNLTE(oo) = 0.0;
  end

