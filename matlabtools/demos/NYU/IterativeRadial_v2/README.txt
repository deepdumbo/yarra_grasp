This package contains the algorithm described in the following publication for iterative reconstruction from undersampled radial MRI using a Total-Variation (TV) constraint:
 
Block KT, Uecker M, Frahm J. Undersampled radial MRI with multiple coils. Iterative image reconstruction using a total variation constraint. Magn Reson Med. 2007 Jun;57(6):1086-98

There are two different Matlab programs: iterativeradial_phantom.m demonstrates the reconstruction of a numerical Shepp-Logan phantom using iterative reconstruction with a TV constraint on the real-part of the image (corresponding to Figure 3 in the publication). iterativeradial_multicoil.m demonstrates the full algorithm with initial coil-estimation step and final joint-coil image calculation (according to Figures 6-8). Real MRI datasets are provided from a phantom and a brain scan, both acquired using the Golden-Angle ordering scheme, which allows retrospective undersampling by truncating spokes at the beginning of the data (note that it is better to truncate at the beginning due to steady-state effects). The optimization is done using the L-BFGS algorithm. The cost functions are implemented in separate files (costfunction_phanom.m for the Shepp-Logan case, and costfunction_coils.m and costfunction_image.m for the full algorithm). The estimate vector x and to-be-calculated gradient vector g are passed as 1D vectors and are reshaped using the helper functions img_to_vec and vec_to_img. In all cost functions, preconditioning with the DCF is used to accelerate convergence. The penalty terms are implemented in external helper functions (TV.m for the Total Variation of the real part, L2D.m for the L2-based smoothness penalty for the coil profiles). These methods have been implemented for simplicity and are not optimized at all for performance. Reconstruction parameters for both Matlab programs are defined at the beginning of both Matlab programs using the structure param. Retrospective undersampling can be defined using the variable param.undersampleSpokesTo. To increase calculation speed, the radial data is downsampled from 2x oversampling prior to the reconstruction if param.readoutDownsampling is set to 1. At the end of the calculation, the iterative reconstruction will be shown in comparison to a gridding solution. Also the estimated coil profiles will be shown. The Shepp-Logan phantom program additionally shows the Fourier transform of both reconstructions, which allows verifying that applying the TV constraint in image space leads to compensation of data gaps in k-space. If the source code is used for reconstruction of other datasets than the provided files (the Matlab program allows reading files from Siemens MR systems with software version VBxx / VDxx), it is necessary to adjust the penalty weights (lambdaxxx) and stopping criteria for the optimizer (stopTolxxx).

Additional comments on the implementation can be found in the source code as well as the original publication. If you are using the code for your research, please cite the publication listed above.

The source code uses the following external packages:
    - NUFFT toolkit by Jeffrey Fessler 
      (http://www.eecs.umich.edu/~fessler/)
    - NUFFT operator by Miki Lustig 
      (http://www.eecs.berkeley.edu/~mlustig/Software.html)
    - Siemens TWIX file reader by Philipp Ehses
    - Poblano Toolbox by Sandia National Laboratories 
      (http://software.sandia.gov/trac/poblano)
    - MRI Phantom by Ronald Ouwekerk
      (http://www.mathworks.com/matlabcentral/fileexchange/1759-mriphantom)

Version 12.10.15. Please check http://cai2r.net/resources/software for updates.
