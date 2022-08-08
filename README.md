# Current status of FastFInS

This is the public version of FastFInS

% run fastfins_check_solvers() to check all the interfaces are installed
* now fem_hp is formally supported by fastfins
* distmesh (by Per-Olof Persson) is used by the spectral FEM package under fem_hp for complex geometries. 
* distmesh is under the terms of the GNU General Public License (version 2 or later). 
* tomo2d contains a function of the tomography code written by J. Heikkinen for computing the cut length. 

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

