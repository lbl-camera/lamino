# Candidate journals — XMCD laminography methods paper

## Paper scope (settled)

A **methods** paper. Contributions:

1. Recoverability characterization of the XMCD forward operator — the
   longitudinal ray transform's null space is the curl-free part of **m**
   (solenoidal injectivity), specialized to the laminographic geometry and to
   minimal azimuthal sampling (two γ orientations).
2. Coupled vector-valued TV via split Bregman on the full 3×3 Jacobian of the
   magnetization field (single isotropic shrinkage over all nine entries,
   `src/bregman.cpp:98-102`).
3. A constant-magnitude constraint |m| = M_s that recovers the null-space
   (curl-free) component a linear inversion provably cannot.
4. Supporting: NUFFT/Toeplitz operator implementation, multi-GPU scaling,
   fiducial bundle-adjustment alignment (`alignlsq`).

**No physical interpretation of any real specimen.** EMMA-665 appears only as a
demonstration that the pipeline runs on measured projections; the science paper
is written separately by the collaborator and will cite this one.

---

## Tier 1 — DOE-countable as high-impact

Impact factors are approximate and move year to year; verify before relying on
them.

### npj Computational Materials — **top pick**
- IF ≈ 10. Nature-family, which is what DOE reporting notices.
- Explicitly a computational-methods venue, so "no science case" is **not** a
  strike against it — unlike the Nature-family journals that require a materials
  result.
- Micromagnetic ground truth plus physics-constrained inversion is dead center
  of scope.
- APC covered by the LBNL/UC Springer Nature read-and-publish agreement.
- Lead with contribution 3 (the constraint defeating a provable null space);
  contributions 1-2 support it.

### Optica
- IF ≈ 9-10. Publishes imaging-methodology firsts with a theoretical spine.
- APC covered by the LBNL Optica agreement.
- Frame as minimal-azimuthal-sampling vector reconstruction. Slightly higher
  risk than npj: expects broad optics interest, and the audience is less
  invested in magnetization specifically.

### Science Advances — only if the constraint result is visually striking
- DOE-countable, broad readership.
- Realistically needs the |m| = M_s recovery to work convincingly on real data,
  not just mumax. Hold in reserve; do not build the manuscript toward it.

---

## Tier 2 — strong, respected, not DOE high-impact

### IEEE Transactions on Computational Imaging
- IF ≈ 4. Best fit if the solver stack is foregrounded: Toeplitz-embedded normal
  operator, split Bregman vs. NAG+qGGMRF, multi-GPU scaling.
- Natural home for the qGGMRF lineage (Bouman–Sauer).
- Will want timing and scaling numbers the current draft lacks.

### Optics Express
- IF ≈ 3-4. Fast turnaround, broad indexing, Optica agreement covers APC.
- Sensible middle path if Tier 1 rejects and speed matters.

---

## Tier 3 — safe landing

### Journal of Synchrotron Radiation
- IF ≈ 2-3. The laminography/beamline community's own venue.
- Marcus, *Dichroic laminography: theory and practice* (JSR 2022) is already in
  the bibliography — this writes directly into that conversation.
- Minimal rework from any Tier 1/2 submission. Reliable acceptance.

### Journal of Applied Crystallography
- Adjacent IUCr option; weaker fit than JSR for magnetic imaging.

---

## Recommended sequence

1. **npj Computational Materials**
2. **Optica**
3. **IEEE TCI**
4. **Journal of Synchrotron Radiation**

---

## Practical notes

- **Post to arXiv early**, before any journal decision. Methods review can run
  6+ months, and the collaborator's science paper must be able to cite this
  regardless of timeline. An arXiv DOI is citable immediately and timestamps the
  methodology.
- **Agree in advance which reconstructions appear in which paper**, so the same
  figure does not surface in both.
- **Zenodo DOIs** for the exact `lamino`, `alignlsq`, and `helmholtz` commits
  used. Helps at review and satisfies the DOE Public Access Plan (PAGES).
- **Phrase the split-Bregman novelty precisely.** 3D TV via split Bregman is
  well published (CT, MRI). The defensible claim is the *coupled 3×3 Jacobian of
  a vector field*, not "split Bregman in 3D" — that phrasing invites a one-line
  desk rejection.
- **Cite the classical null-space result** (Sharafutdinov 1994;
  Phatak/Beleggia/De Graef on vector-field electron tomography) rather than
  presenting the div-free/curl-free split as new. The novelty is defeating the
  null space, not finding it. Verify both references before citing — they are
  from memory.
