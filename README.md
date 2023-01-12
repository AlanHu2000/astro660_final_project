# astro660_final_project
Reproduce dust continuum map and polarization map of HL tau disk through RADMC3D.

## Step 0 : Install RADMC3D
If your machine does not install RADMC3D, please download and install at : \
https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/

## Step 1 : Run Self-Scattering produced polarization image
**1.Change the problem setting for RADMC3D** \
Edit "./dust_scat_pol_image/problem_setup_dust_scat.py" to change disk parameters, grain size or running setting for RADMC3D 

**2.Generate the input files for RADMC3D** \
Execute "./dust_scat_pol_image/problem_setup_dust_scat.py" to generate input files for RADMC3D. 

**3.Run RADMC3D calculation and generate output image** \
Type "radmc3d image lambda 870 incl 46.7 npix 300 posang 42 dpc 140 stokes" to generate full stokes image of HL Tau at ALMA band 7 $(870\mu m)$. \
Type "radmc3d tausurf 0.1 lambda 870 incl 46.7 npix 300 dpc 140" to generate optical depth $\tau=0.1$ surface.

**4.Optional**\
If you want to run the Self-Scattering of different maximum grain size, you can change file "./dust_scat_pol_imagedustkapscatmat_new_100" to any files in "./scatmat_max_birnstiel2018/", do not forget to change "./dust_scat_pol_image/problem_setup_dust_scat.py" as well.

## Step 2 : Run Grain Alignment produced polarization image
**1.Run RADMC3D** \
Change "./dust_scat_pol_image/problem_setup_dust_scat.py" to "./grain_align_pol_image/problem_setup_alignment_oblate.py" and do the Step 1-1 to 1-3 again.

**2.Optional** \
If you want to try different axis ratio of oblate grains or alignment efficiency, please edit "./grain_align_pol_image/problem_setup_alignment_oblate.py" to change "./grain_align_pol_image/grainalign_dir.inp".

## Step 3 : Visualization of polarization image
Use jupyter notebook "draw_image.ipynb" file to do the visualization.
