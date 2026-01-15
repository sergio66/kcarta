%% user has set his/her inputs
user_set_input_opticaldepths

%% user has set paths to files
user_set_dirs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

cd ..

iDownLook = +1;
initialize_extra;
[h,ha,p,pa] = doload(dirin,fin,klayers_code,iAirs); %%% load in LAYERS profile

global iDebug
%%% +10XY for Tz debug, 20XY for Qz debug, 1 for talky, 0 quiet as a mouse
%iDebug = +1;      
%iDebug = +1090;   
iDebug = +0;      

for ip = iProfRun
  fprintf(1,'processing profile %5i \n',ip);
  [head, prof] = subset_rtp(h, p, h.glist, [], ip);
  [nlays,prof,rFracBot,ropt0] = initialize_kcmix(prof,iDownLook,ropt0);
  aux_struct = auxiliary_set(fA,fB,nlays,rFracBot,CKD,cswt,cfwt,refp);
  odsOut = opticaldepths(head,prof,aux_struct,ropt0,iGasDoOD,iBreakoutCont);
end

stuff.freqs = [fA fB]; stuff.input_rtpfile = [dirin '/' fin];
stuff.layersprof = prof;  stuff.iGasDoOD = iGasDoOD;   stuff.iProfRun = iProfRun;
stuff.iBreakoutCont = iBreakoutCont;

clear aux_struct iProfRun iHITRAN
clear head prof ha pa CKD cdir co2ChiFilePath iAirs iDebug iDownLook 
clear ip ropt refp cswt cfwt
clear junkdir klayers_code kpath* nlays nltedir rFracBot soldir
clear fA fB dirin fin kdatadir iMatlab_vs_f77
clear h p iDoJac iJacobOutput  str str0 iBreakoutCont iGasDoOD
clear iBreakoutCont iGasDoOD

cd Test

toc
