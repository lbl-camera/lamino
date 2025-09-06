
#include "array.h"
#include "nufft.h"

#include "optimize.h"

template <typename Real_t>
Array<Real_t> run_mbir(Array<Real_t> sinogram, const vector<Real_t> &theta) {

    int nprojs = sinogram.nrows();
    int ncols = sinogram.ncols();

    // inital setup
    auto bT = backproject(sinogram, theta);

    // form the polar grid
    auto pg = nufft::polar_grid<Real_t>(theta, sinogram.ncols());

    // calculate psf
    auto psf = calc_psf(pg);

    // create callable objects for the optimization
    auto cost_function = [&](const Array<Real_t> &x) {
        return (forward(x, theta) - sinogram).norm();
        };
    auto cost_gradient = [&](const Array<Real_t> &x) {
        return gradient(x, pg) - bT;
        };

    // create optimization object
    auto o = opt::Optimizer<Real_t, Array, decltype(cost_function), decltype(cost_gradient)>
        (cost_function, cost_gradient);

    // run the optimization
    int num_iters = 100;
    Real_t tol = 1.0e-04;
    Real_t xtol = 1.0e-03;

    // Lipschitz constant
    Array<Real_t> x0(ncols, ncols, 1);
    auto L = (forward(x0, theta)).max();

    return o.run(x0, num_iters, L, tol, xtol);
}
