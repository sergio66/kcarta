addpath USER_SET

%% user has set his/her inputs
%user_set_input_breno
%user_set_input_downlook

user_set_input_downlook_testjacs_Jan2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% user has set paths to files
%user_set_dirs
user_set_dirs_Jan2026
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
%iDebug = +1;      

for ip = iProfRun
  fprintf(1,'processing profile %5i \n',ip);
  [head, prof] = subset_rtp(h, p, h.glist, [], ip);
  [nlays,prof,rFracBot,ropt0] = initialize_kcmix(prof,iDownLook,ropt0);
  aux_struct = auxiliary_set(fA,fB,nlays,rFracBot,CKD,cswt,cfwt,refp);
  if exist('gasidlist')
    disp('gasidlist exists ... so using that')
    aux_struct.gasidlist = gasidlist;
  end
  if iDoJac(1) > 0
    [radsOut,jacsOut,gasinfo_per_chunk] = downlook_jac(head,prof,aux_struct,ropt0,iDoJac,iJacobOutput);
  else
    [radsOut] = downlook_nojac(head,prof,aux_struct,ropt0);
    end
  end

if iDoJac(1) > 0
  stuff.gasinfo_per_chunk = gasinfo_per_chunk;
  %% so you can look at eg stuff.gasinfo_per_chunk.info{4} and it will tell you wchunk followed by gases used
end
stuff.freqs = [fA fB]; stuff.input_rtpfile = [dirin '/' fin];
stuff.layersprof = prof;  stuff.iDoJac = iDoJac;   stuff.iProfRun = iProfRun;
stuff.iJacobOutput = iJacobOutput;

clear aux_struct iProfRun iHITRAN
clear head prof ha pa CKD cdir co2ChiFilePath iAirs iDebug iDownLook 
clear ip ropt refp cswt cfwt
clear junkdir klayers_code kpath* nlays nltedir rFracBot soldir
clear fA fB dirin fin kdatadir iMatlab_vs_f77
clear h p iDoJac iJacobOutput str str0

cd Test

toc
