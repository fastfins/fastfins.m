# Current status of FastFInS

% run fastfins_check_solvers() to check all the interfaces are installed
* now fem_hp is formally supported by fastfins
* distmesh is used by the spectral FEM package under fem_hp, need to check the license for distributing this package
* tomo2d uses the tomography code written by J. Heikkinen to setup geometry of the forward model. We have the permission to use it. This may need to be re-written before releasing the code.

# Examples

The following examples need to added (under the examples)

* DILI and H-MALA
* RTO
* Poisson with tensor field
* p-Poisson 
* Stokes equation (using the spectral code)

More work needs:

* RTO needs to be updated (explicitJ and trust region)
* optimisation solver needs more work: newton CG and trust region search

The following package needs to be removed (or significantly updated)
* fem_bilinear
