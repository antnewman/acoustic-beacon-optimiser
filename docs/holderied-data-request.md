# Draft email: Simon et al. (2020) FEM data request

This is a draft email to send to the senior author of Simon et al. (2020,
PNAS 117: 1367-1374) requesting the raw FEM scattering results, which are
needed for a cross-validation step in the `acoustic-beacon-optimiser`
paper. Review and edit before sending.

**To:** `marc.holderied@bristol.ac.uk` (publicly listed contact; confirm
via Bristol staff directory before sending)

**Cc:** `ralph.simon@fu-berlin.de` (first author; include only if you want
Ralph Simon copied in)

**Subject:** `FEM data request: Simon et al. 2020 PNAS for open-source BEM validation`

---

## Draft body

Dear Professor Holderied,

I am writing to ask whether it would be possible to share the raw
finite-element scattering results from your 2020 PNAS paper
(Simon et al., "Biosonar-inspired acoustic reflectors increase detection
distances of small objects", DOI: 10.1073/pnas.1909890117) with a view to
cross-validating an independent numerical solver.

Briefly, I am developing an open-source Python package,
`acoustic-beacon-optimiser`, that applies the Helmholtz boundary element
method (Burton-Miller CFIE via `bempp-cl`) to floral acoustic reflectors
and performs formal shape optimisation against biologically meaningful
objective functions. To my knowledge this is the first BEM computation
for this class of reflector and the first systematic shape-optimisation
study of bat-pollinated floral geometries.

The solver has been validated against the analytical Mie series for a
rigid sphere, with agreement to approximately 0.03 dB at Helmholtz
numbers He = 5 and He = 10. For a paper-grade numerical validation on
a relevant geometry I would like to reproduce the nine spherical-cap
configurations in your 2020 paper (r = 35, 50, 70 mm at d = 25, 30,
49 mm, where geometrically valid) and compare BEM target strength
against the FEM values you computed.

If you are willing to share any of the following I would be very
grateful:

1. Target-strength tables or spectra (frequency x incidence angle) for
   the nine spherical-cap configurations.
2. Impulse responses (or equivalent) for those configurations, if
   different from the above.
3. Any mesh files or geometry definitions you used, for a like-for-like
   numerical comparison.

I would of course cite the data appropriately in any resulting
publication, and I am happy to send the draft manuscript and validation
notebook for your comments before submission.

The repository (MIT licence, early development) is at:
[https://github.com/antnewman/acoustic-beacon-optimiser](https://github.com/antnewman/acoustic-beacon-optimiser)

Please let me know if there is any institutional or ethical process I
should follow to make this request formally.

With many thanks and kind regards,

Ant Newman, CEng MIET MIEEE
[tortoiseai.co.uk/about/ant-newman](https://tortoiseai.co.uk/about/ant-newman)
`antjsnewman@outlook.com`

---

## Notes before sending

- Confirm Marc Holderied's current email via the Bristol staff directory
  rather than trusting the `@bristol.ac.uk` guess.
- Decide whether to copy Ralph Simon: he was first author in 2011
  (*Science*) and co-author in 2020, and runs a related research line
  at FU Berlin; copying him increases chances of a reply but also
  makes the thread less private.
- If you want to offer authorship or formal acknowledgement in exchange
  for raw data, add a short paragraph on that. Current draft offers
  citation only, which is the lightest-touch ask.
- Avoid making implicit claims about when the paper will be submitted
  or accepted. The draft above stays non-committal on timelines.
- Do not attach the draft manuscript on first contact; offer it.
  Sending large PDFs unsolicited is ungraceful.
