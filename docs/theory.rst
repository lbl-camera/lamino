Theory
======

A forward and backward projection operators are fundamental components of the iterative reconstruction process. tomoCAM-lamino implements these operators efficiently using non-uniform fast Fourier transforms (NUFFT). 

.. math::
    \mathcal(P) \vec{m} = \int_{-\infty}^\infty \vec{m} \cdot \hat{e} \, e^{2 \pi i \vec{q} \cdot r \hat{e}} \, dr