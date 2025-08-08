%-- 6/18/2025 12:59 AM --%
cd('/Volumes/PrecunealSSD2/AnalyticSignalHCP');
imX1 = cifti_read("imX_as_sub-n8_ses-n8_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-fft-ASHCPPar-20250618002801_par1.dtseries.nii")
imX2 = cifti_read("imX_as_sub-n8_ses-n8_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-fft-ASHCPPar-20250618002733_par2.dtseries.nii")
imagesc(imX1.cdata'); colorbar; clim([-200, 200]);
figure; imagesc(imX2.cdata'); colorbar; clim([-200, 200]);
clim([-50, 50])
%-- 6/20/2025 10:56 PM --%
as = mlraut.AnalyticSignalHCPPar.load("sub-993675_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z = as.zeta(as.bold_signal, physio_signal);
z = as.zeta(as.bold_signal, as.physio_signal);
open mlraut.Cifti
as.cifti.write_cifti(real(z), "rezeta_as_sub-993675_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as.test_ref_dscalar_fqfn
as.cifti.test_ref_dscalar_fqfn
as.cifti.write_cifti(real(z), "rezeta_as_sub-993675_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
mg
this.task_dir
as = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z = as.zeta(as.bold_signal, as.physio_signal);
imagesc(real(z)); colorbar
clim([-3,3])
figure; imagesc(imag(z)); colorbar; clim([-3,3])
as.cifti.write_cifti(real(z), "rezeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as.cifti.write_cifti(imag(z), "imzeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
as = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-brightest-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z = as.zeta(as.bold_signal, as.physio_signal);
figure; imagesc(real(z)); colorbar; clim([-3,3])
figure; imagesc(imag(z)); colorbar; clim([-3,3])
figure; imagesc(real(as.bold)); colorbar; clim([-3,3])
figure; imagesc(real(as.bold_signal)); colorbar; clim([-3,3])
as_gsr0 = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z_gsr1 = as_gsr1.zeta(as_gsr1.bold_signal, hilbert(as_gsr1.physio_supplementary("iFV")));
z_gsr0 = as_gsr0.zeta(as_gsr0.bold_signal, as_gsr0.physio_signal);
figure; imagesc(real(z_gsr1)); colorbar; clim([-3,3])
figure; imagesc(imag(z_gsr1)); colorbar; clim([-3,3])
figure; imagesc(real(as_gsr0.bold_signal)); colorbar; clim([-3,3])
figure; imagesc(real(as_gsr1.bold_signal)); colorbar; clim([-3,3])
as.cifti.write_cifti(real(z), "rezeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as_gsr1.cifti.write_cifti(real(z_gsr1), "rezeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as_gsr1.cifti.write_cifti(imag(z_gsr1), "imzeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
isnumeric([])
as_gsr1 = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-brightest-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z_gsr1 = as_gsr1.zeta(as_gsr1.neg_dbold_dt(as_gsr1.bold_signal), hilbert(as_gsr1.physio_supplementary("iFV")));
z_gsr1 = as_gsr1.zeta(as_gsr1.neg_dbold_dt(as_gsr1.bold_signal, hilbert(as_gsr1.physio_supplementary("iFV"))), hilbert(as_gsr1.physio_supplementary("iFV")));
iFV_gsr1 = hilbert(as_gsr1.physio_supplementary("iFV")); z_gsr1 = as_gsr1.zeta(as_gsr1.neg_dbold_dt(as_gsr1.bold_signal, iFV_gsr1, iFV_gsr1);
iFV_gsr1 = hilbert(as_gsr1.physio_supplementary("iFV")); z_gsr1 = as_gsr1.zeta(as_gsr1.neg_dbold_dt(as_gsr1.bold_signal, iFV_gsr1, iFV_gsr1));
z_gsr1 = as_gsr1.zeta(as_gsr1.neg_dbold_dt(as_gsr1.bold_signal, iFV_gsr1), iFV_gsr1);
as_gsr1 = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-brightest-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z_gsr1 = as_gsr1.zeta(as_gsr1.bold_signal, iFV_gsr1);
as_gsr1.cifti.write_cifti(real(z_gsr1), "rezeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as_gsr1.cifti.write_cifti(imag(z_gsr1), "imzeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as_gsr0 = mlraut.AnalyticSignalHCPPar.load("sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
z_gsr0 = as_gsr0.zeta(as_gsr0.bold_signal, as_gsr0.physio_signal);
as_gsr0.cifti.write_cifti(real(z_gsr0), "rezeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
as_gsr0.cifti.write_cifti(imag(z_gsr0), "imzeta_as_sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.dtseries.nii")
figure; imagesc(real(z_gsr1)); colorbar; clim([-3,3])
figure; imagesc(real(as_gsr1.neg_dbold_dt(as_gsr1.bold_signal))); colorbar; clim([-3,3])
figure; imagesc(real(z_gsr0)); colorbar; clim([-3,3])
figure; imagesc(real(as_gsr0.neg_dbold_dt(as_gsr0.bold_signal))); colorbar; clim([-3,3])
ls = load("fftzeta_as_sub-n4_ses-n4_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20250621160915.mat")
ld
figure; imagesc(real(ld.zeta_)); colorbar
clim([-100,100])
clim([-30,30])
clim([-300,300])
ld = load("fftzeta_as_sub-n4_ses-n4_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20250621160915.mat")
figure; imagesc(real(ld.zeta_)); colorbar
clim([-300,300])
clim([-30,30])
clim([-10,10])
ld = load("fftzeta_as_sub-n4_ses-n4_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20250621160915.mat")
figure; imagesc(real(ld.zeta_)); colorbar
clim([-10,10])
ld = load("fftzeta_as_sub-n4_ses-n4_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20250621160915.mat")
figure; imagesc(real(ld.zeta_)); colorbar
clim([-10,10])
clim([-10,10])
clim([-1000,1000])
clim([-100,100])
clim([-30,30])
ld = load("fftzeta_as_sub-n4_ses-n4_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar-20250621170613.mat")
figure; imagesc(real(ld.zeta_)); colorbar
clim([-1e2,1e2])
clim([-1e3,1e3])
clim([-500,500])
clim([-1e3,1e3])
figure; imagesc(real(ifft(ld.zeta_))); colorbar
clim([-100, 100])
clim([-30, 30])
globbed = mglob("9*/fftzeta_as_sub-n4_*.mat");
for g = globbed
ld = load(g);
for g = globbed(1:5:end)
ld = load(g);
figure; imagesc(real(ld.zeta_)); colorbar; end
clim([-100,100])
clim([-300,300])
clim([-1000,1000])
clim([-300,300])
for g = globbed(1:5:end)
ld = load(g);
figure; imagesc(real(ld.zeta_)); colorbar; clim([-500,500]); end
clim([-100,100])
clim([-1000,1000])
for idx = 1:5
for g = globbed((idx-1)*10+1:idx*10)
ld = load(g);
figure; imagesc(real(ld.zeta_)); colorbar; clim([-500,500]); end
for g = globbed
ld = load(g);
figure; imagesc(real(ld.zeta_)); colorbar; clim([-500,500]); saveFigure2(gcf, mybasename(g), closeFigure=true); end
for g = globbed
fp = strrep(g, filesep, "_");
xlabel("greyordinate")
figure; imagesc(real(fftshift(ld.zeta_, 1)));
clim([-500,500])
as = mlraut.AnalyticSignalHCPPar.load("996782/sub-996782_ses-rfMRI-REST1-LR_proc-iFV-brightest-gsr1-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
as = mlraut.AnalyticSignalHCPPar.load("996782/sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
as.rugplot_fft_zeta(ld.zeta_)
as = mlraut.AnalyticSignalHCPPar.load("996782/sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
as.rugplot_fft_zeta(ld.zeta_)
as = mlraut.AnalyticSignalHCPPar.load("996782/sub-996782_ses-rfMRI-REST1-LR_proc-iFV-gsr0-ddt1-butter8-lp0p1-hp0p01-scaleiqr-subset-ASHCPPar.mat", class="mlraut.AnalyticSignalHCPPar")
as.rugplot_fft_zeta(ld.zeta_, fig_fileprefix="test")