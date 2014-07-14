
% compare 100 layer single gas absorptions from 
% kcmix and the new kcarta (v 1.07)
%
% kcarta data is from functions "onegas" and "xsecgas"
% in the directory kctest

% path to kcarta scripts
addpath /home/motteler/radtrans/kctest

% chunk start wavenumber
vchunk = input('chunk wavenumber > ');

% select a point profile (needs a full path)
% ptpro = input('GENLN2 format point profile > ', 's');
% ptpro = '/home/motteler/radtrans/klayers/Data/Pexample/myp1';
ptpro = '/home/motteler/kcmix2/refpro.ptpro';

% read point profile as a kcmix structure
p1 = pt2kcmix(ptpro);

% specify gas set
gid = input('HITRAN gas ID > ');
gid = intersect(gid, p1.glist);
gind = find(p1.glist == gid);

% select specified gas
p1.glist = p1.glist(gind);
p1.gamnt = p1.gamnt(:,gind);
p1.gpart = p1.gpart(:,gind);

% location of compressed data 
cdir = '/asl/data/kcarta/v20.matlab';

% calculate new absorptions with kcmix
fprintf(1, 'calling kcmix...\n');
[a1, freq] = kcmix(p1, vchunk, cdir);

% compare to kcarta tabulation
fprintf(1, 'calling kcarta...\n');
if gid <= 50
  [kabsc, kfreq] = onegas(ptpro, vchunk, gid);
else
  [kabsc, kfreq] = xsecgas(ptpro, vchunk, gid);
end

% compare total column absorption

% figure(2); clf
% subplot(2,1,1)
% plot(freq, sum(a1'), kfreq, sum(kabsc'));
% legend('100 layer kcmix', '100 layer kcarta')
% title('total column absorption')
% 
% subplot(2,1,2)
% plot(freq, (sum(a1') - sum(kabsc')));
% title('difference')
% % plot(freq, (sum(a1') - sum(kabsc')) ./ sum(a1') );
% % title('relative difference')

% compare selected layer absorptions

tlay = 40; 
while  ~isempty(tlay)

  figure(1); clf
  subplot(2,1,1)
  plot(freq, a1(:,tlay), kfreq, kabsc(:,tlay));
  title(sprintf('gas %d absorption at layer %d', gid, tlay))
  legend('100 layer kcmix', '100 layer kcarta')

  subplot(2,1,2)
  % plot(freq, (a1(:,tlay) - kabsc(:,tlay)));
  % title('difference')
  plot(freq, (a1(:,tlay) - kabsc(:,tlay)) ./ a1(:,tlay));
  title('relative difference')

  tlay = input('test layer > ');
end
