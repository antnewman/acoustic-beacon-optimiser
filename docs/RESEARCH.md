# Mathematics of bat-flower acoustic coevolution: literature review and gap analysis

> Last updated: 2026-04-14

This document surveys the quantitative literature on acoustic interactions between
chiropterophilous flowers and their glossophagine pollinators, identifies the
mathematical tools already applied, and catalogues seven tractable open problems
where new modelling contributions would advance the field.

---

## 1 Floral acoustic reflectors

Four distinct reflector morphologies have been documented, each exploiting a
different physical mechanism. Together they span the Helmholtz number range
He = ka ~ 3-30 (where k is the acoustic wavenumber and a a characteristic
flower dimension at 50-130 kHz), placing them squarely in the resonance
scattering regime where neither geometric optics nor Rayleigh scattering
applies. Despite this, no study has performed shape optimisation over any
biologically meaningful objective function.

### 1.1 Mucuna holtonii vexillum (cat's-eye retroreflector)

The upright petal (vexillum) of *Mucuna holtonii* was the first floral
structure shown to function as an acoustic guide for bat pollinators.
von Helversen and von Helversen (1999) demonstrated that the concave vexillum
acts as a cat's-eye retroreflector, returning a conspicuous echo to approaching
bats over a wide angular range. Echo discrimination operates at distances of
15-50 cm. In a follow-up study, von Helversen and von Helversen (2003)
characterised the spectral composition of reflected echoes and showed that it
varies systematically with incidence angle, providing bats with both range and
angular cues.

Critically, no scattering model (boundary element, physical optics, or
otherwise) has been fitted to the vexillum geometry. The published analyses
remain purely empirical.

- von Helversen, D. and von Helversen, O. (1999). Acoustic guide in
  bat-pollinated flower. *Nature*, 398, 759-760.
  DOI: [10.1038/19648](https://doi.org/10.1038/19648)
- von Helversen, D. and von Helversen, O. (2003). Object recognition by
  echolocation: a nectar-feeding bat exploiting the flowers of a
  rain forest vine. *Journal of Comparative Physiology A*, 189, 327-336.
  DOI: [10.1007/s00359-003-0405-3](https://doi.org/10.1007/s00359-003-0405-3)

### 1.2 Marcgravia evenia dish-leaf (multidirectional beacon)

The concave leaf of *Marcgravia evenia* represents the most thoroughly
characterised floral acoustic reflector. Simon et al. (2011) showed that the
dish-shaped leaf above the inflorescence acts as a multidirectional acoustic
beacon, halving the foraging time of *Glossophaga soricina* in flight-cage
experiments. Ensonification was performed using maximum-length sequences (MLS)
with impulse response (IR) extraction.

Simon et al. (2020) extended this work with a finite-element method (FEM)
simulation of spherical caps, providing the only published numerical model of a
floral reflector. Key results include a target strength (TS) formula relating
cap geometry to echo amplitude:

> TS = 20 log10(a^2 / (4 lambda d))

where a is the cap aperture radius, lambda the wavelength, and d the focal
distance. The impulse response exhibits a two-peak structure whose temporal
separation encodes reflector depth via the path-length difference between the
specular return and the rim-diffracted wave.

- Simon, R., Holderied, M. W., Koch, C. U. and von Helversen, O. (2011).
  Floral acoustics: conspicuous echoes of a dish-shaped leaf attract bat
  pollinators. *Science*, 333, 631-633.
  DOI: [10.1126/science.1204210](https://doi.org/10.1126/science.1204210)
- Simon, R., Rupitsch, S. J. and Baumann, M. J. (2020). Numerical simulation
  of bat-pollinated flower reflectors. *Proceedings of the National Academy of
  Sciences*, 117, 1367-1374.
  DOI: [10.1073/pnas.1909890117](https://doi.org/10.1073/pnas.1909890117)

### 1.3 Espostoa frutescens cephalium (absorption contrast)

Rather than enhancing echoes, the woolly cephalium of the cactus *Espostoa
frutescens* operates by a contrasting mechanism: it absorbs ultrasound, creating
a "dark" acoustic marker against the bright specular background of the smooth
cactus body. Simon et al. (2023) measured -14 dB attenuation at 90 kHz, a
frequency closely matching the dominant harmonic of the pollinating bat's
echolocation call.

No porous media acoustic model has been applied to the cephalium. The Delany-
Bazley empirical model for fibrous absorbers or, for a more physically grounded
approach, Biot poroelastic theory would be appropriate starting points. The
trichome density and fibre diameter data required for parameterisation are
available from scanning electron micrographs in the original publication.

- Simon, R., Knornschild, M., Tschapka, M., Grund, A., Kalko, E. K. V. and
  Hoelzle, L. E. (2023). Acoustic absorption by cactus trichomes facilitates
  bat pollination. *Journal of Experimental Biology*, 226, jeb245263.
  DOI: [10.1242/jeb.245263](https://doi.org/10.1242/jeb.245263)

### 1.4 Bell-shaped flowers (spectral directivity)

von Helversen, Holderied and von Helversen (2003) measured the echo spectra of
bell-shaped flowers from multiple incidence angles and characterised the
resulting spectral directional patterns. These patterns are, in physical
acoustics terminology, bistatic scattering cross-sections parameterised by
frequency and angle.

The analysis was entirely empirical. No analytical model (e.g., cylindrical
shell scattering, T-matrix formulation) was fitted to the measured directivity
patterns, despite the relatively simple axisymmetric geometry of bell-shaped
corollas.

- von Helversen, D., Holderied, M. W. and von Helversen, O. (2003). Echoes
  of bat-pollinated bell-shaped flowers: conspicuous for nectar-feeding bats?
  *Journal of Experimental Biology*, 206, 1025-1034.
  DOI: [10.1242/jeb.00203](https://doi.org/10.1242/jeb.00203)

---

## 2 Glossophagine echolocation signal design

### 2.1 Search-approach-terminal phase structure

Flower-visiting bats modulate their echolocation calls in a stereotyped sequence
as they approach a target, analogous to the well-known search-approach-terminal
(SAT) phases of insectivorous bats. Gonzalez-Terrazas et al. (2016) documented
this sequence in *Leptonycteris yerbabuenae* approaching columnar cacti. Pulse
duration shortens, inter-pulse interval (IPI) decreases, and bandwidth increases
as the bat transitions from search to terminal phase. IPI sets the maximum
unambiguous range via the relationship R_max = c * IPI / 2; the systematic
decrease in IPI therefore trades range coverage for update rate as the bat
closes on the flower.

- Gonzalez-Terrazas, T. P., Koblitz, J. C., Fleming, T. H., Medellein, R. A.,
  Kalko, E. K. V., Schnitzler, H.-U. and Tschapka, M. (2016). How
  nectar-feeding bats localize their food: echolocation behaviour and spatial
  memory in *Leptonycteris yerbabuenae*. *PLoS ONE*, 11, e0163492.
  DOI: [10.1371/journal.pone.0163492](https://doi.org/10.1371/journal.pone.0163492)

### 2.2 Glossophaga soricina call parameters

Simon et al. (2006) characterised the echolocation calls of *Glossophaga
soricina* during flower approach. Calls are frequency-modulated (FM) sweeps
spanning two harmonics from 56 to 137 kHz, with notably low source levels
compared to insectivorous species. The low intensity is consistent with the
short operating range (typically < 1 m) and the need to avoid alerting competing
bats.

- Simon, R., Holderied, M. W. and von Helversen, O. (2006). Size
  discrimination of hollow hemispheres by echolocation in a nectar-feeding bat.
  *Journal of Experimental Biology*, 209, 3599-3609.
  DOI: [10.1242/jeb.02398](https://doi.org/10.1242/jeb.02398)

### 2.3 FM sweep instantaneous frequency models

The instantaneous frequency f(t) of a bat FM sweep has been modelled using
several functional forms:

- **Hyperbolic**: f(t) = k / (t + t_0). Altes and Titlebaum (1970) showed that
  hyperbolic sweeps are Doppler-tolerant, meaning the matched-filter output is
  insensitive to target radial velocity. This is the theoretically optimal
  waveform for simultaneous range-Doppler estimation.
  DOI: [10.1121/1.1912222](https://doi.org/10.1121/1.1912222)

- **Exponential**: f(t) = f_0 * exp(-alpha * t). Boonman and Schnitzler (2005)
  fitted exponential models to the downward FM sweeps of several vespertilionid
  species and found good empirical agreement.
  DOI: [10.1007/s00359-004-0566-8](https://doi.org/10.1007/s00359-004-0566-8)

- **Logarithmic**: f(t) = f_0 + beta * ln(t / t_0). Masters, Jacobs and
  Simmons (1991) proposed logarithmic models based on cochlear representations,
  arguing that uniform sampling on a log-frequency axis produces constant neural
  delay increments across the tonotopic map.

No study has parameterised these models specifically for glossophagine FM sweeps,
nor compared the goodness-of-fit across models for this clade. This constitutes
a tractable gap (see Section 6.4).

### 2.4 Adaptive sensing and proportional navigation

Echolocation is not merely passive reception; bats actively control their
emission parameters to optimise information gain. Ghose and Moss (2006a)
demonstrated that *Eptesicus fuscus* uses a proportional navigation strategy
during insect interception, adjusting flight trajectory based on the rate of
change of the bearing angle to the target. Ghose et al. (2006b) further showed
that sonar beam aim, head orientation, and pulse timing are coordinated as part
of an integrated adaptive sensing system.

Although these studies concern insectivorous species, the principles of adaptive
sensing apply directly to flower approach, where the bat must control its
trajectory to dock with a small, stationary nectar source.

- Ghose, K. and Moss, C. F. (2006). Steering by hearing: a bat's acoustic
  gaze is linked to its flight motor output by a delayed, adaptive linear law.
  *Journal of Neuroscience*, 26, 1704-1710.
  DOI: [10.1523/JNEUROSCI.4315-05.2006](https://doi.org/10.1523/JNEUROSCI.4315-05.2006)
- Ghose, K., Horiuchi, T. K., Krishnaprasad, P. S. and Moss, C. F. (2006).
  Echolocating bats use a nearly time-optimal strategy to intercept prey.
  *PLoS Biology*, 4, e108.
  DOI: [10.1371/journal.pbio.0040108](https://doi.org/10.1371/journal.pbio.0040108)

### 2.5 Neural network flower-docking model

Nguyen and Vanderelst (2025) trained a neural network to control a simulated
bat's approach to a flower, demonstrating that end-to-end learning can replicate
the docking manoeuvre without explicit modelling of the echo processing pipeline.
While this represents a useful proof of concept, the black-box nature of the
model limits its explanatory power. The approach does not decompose the problem
into perception and control components and therefore cannot isolate the
contribution of floral reflector geometry to docking performance.

### 2.6 Missing formalism: POMDP for glossophagine flower approach

No published study formulates the glossophagine flower approach as a partially
observable Markov decision process (POMDP), despite the natural fit: the bat
maintains a belief state over flower position (partially observable from noisy
echoes), takes actions (flight manoeuvres and call emissions), and receives
stochastic observations (echo returns). A POMDP formulation would provide a
principled framework for understanding adaptive sensing during approach (see
Section 6.6).

---

## 3 Echo processing models

### 3.1 SCAT model

The Spectrogram Correlation and Transformation (SCAT) model is the dominant
computational framework for bat echo processing. Originally proposed by Saillant
et al. (1993), it implements a three-stage pipeline:

1. **Cochlear filterbank**: a bank of gammatone filters simulating the basilar
   membrane, producing a spectrogram-like representation.
2. **Spectrogram correlation**: cross-correlation of the outgoing call
   spectrogram with the returning echo spectrogram to extract delay information.
3. **Spectrogram transformation**: transformation of the correlation output into
   a target-range profile.

Ming et al. (2021) extended the SCAT model and validated it against neural
recordings from the inferior colliculus. Open-source implementations are
available: SCAT (github.com/gomingchen/SCAT) and BiSCAT for bistatic
configurations (github.com/gaudetteje/biscat).

- Saillant, P. A., Simmons, J. A., Dear, S. P. and McMullen, T. A. (1993). A
  computational model of echo processing and acoustic imaging in
  frequency-modulated echolocating bats: the spectrogram correlation and
  transformation receiver. *Journal of the Acoustical Society of America*, 94,
  2691-2712.
  DOI: [10.1121/1.407353](https://doi.org/10.1121/1.407353)
- Ming, C., Haro, S., Bhatt, A. and Bhatt, D. (2021). A neural network model
  for echo processing in the inferior colliculus of echolocating bats. *PLoS
  Computational Biology*, 17, e1008677.
  DOI: [10.1371/journal.pcbi.1008677](https://doi.org/10.1371/journal.pcbi.1008677)

### 3.2 Size-independent shape recognition via log-frequency invariance

A remarkable feature of bat echo classification is its invariance to object
size. Von Helversen (2004) showed that *G. soricina* can discriminate hollow
hemisphere shapes regardless of absolute size, suggesting that the relevant
perceptual features are encoded on a log-frequency axis. When an object is
scaled by factor s, its echo spectrum shifts by log(s) on a log-frequency axis
but retains its shape. Because the mammalian cochlea maps frequency
tonotopically (approximately logarithmically), this spectral shift corresponds
to a spatial translation along the basilar membrane, which the auditory system
can factor out.

Simon et al. (2006) exploited the same principle, demonstrating that size
discrimination of hollow hemispheres depends on spectral shape rather than
spectral position.

- Von Helversen, D. (2004). Object classification by echolocation in nectar
  feeding bats: size-independent generalization of shape. *Journal of
  Comparative Physiology A*, 190, 515-521.
  DOI: [10.1007/s00359-004-0492-9](https://doi.org/10.1007/s00359-004-0492-9)

### 3.3 Machine learning classification of echo spectrograms

Two main machine learning approaches have been applied to bat echo
classification:

**SVM on echo spectrograms.** Yovel et al. (2008) trained support vector
machines on spectrographic representations of echoes from natural objects,
achieving approximately 90% classification accuracy. This demonstrated that
echo spectrograms contain sufficient information for reliable object
discrimination without requiring biologically detailed processing models.

- Yovel, Y., Franz, M. O., Stilz, P. and Schnitzler, H.-U. (2008).
  Plant classification from bat-like echolocation signals. *PLoS Computational
  Biology*, 4, e1000032.
  DOI: [10.1371/journal.pcbi.1000032](https://doi.org/10.1371/journal.pcbi.1000032)

**CNN across multiple species.** Simon et al. (2021) trained convolutional
neural networks on echoes from 12 plant species, achieving 84% accuracy.
The CNN approach captured cross-species variation and demonstrated that echo
classification generalises across the botanical diversity encountered by
foraging bats.

- Simon, R., Rupitsch, S. J., Guenther, F. and Kuc, R. (2021). Echo-based
  plant classification by echolocating bats: a convolutional neural network
  approach. *PLoS Computational Biology*, 17, e1009706.
  DOI: [10.1371/journal.pcbi.1009706](https://doi.org/10.1371/journal.pcbi.1009706)

### 3.4 Missing: wavelet and scattering-transform analysis

Neither wavelet analysis nor scattering transforms (Mallat, 2012) have been
applied to flower echoes. The scattering transform is particularly promising
because it is designed to extract features that are invariant to translation and
stable to deformations; properties closely matching the log-frequency shift
invariance observed in bat echo perception (Section 3.2). This gap is addressed
in Section 6.7.

---

## 4 Coevolutionary modelling: the largest gap

**No published study formally models bat-flower coevolution with equations that
incorporate the acoustic channel.** This is the central gap in the literature.
Pollination network studies treat interaction strengths as abstract quantities;
acoustic studies treat flower geometry as fixed. No framework bridges the two.

Several mathematical tools are available but have not been applied to this
system:

### 4.1 Adaptive dynamics canonical equation

Dieckmann and Law (1996) derived the canonical equation of adaptive dynamics,
which describes the expected rate of phenotypic change in a monomorphic
population under rare mutations:

> dx/dt = (1/2) * mu(x) * n(x) * sigma^2(x) * (partial/partial x') s(x', x) |_{x'=x}

where mu is the mutation rate, n the equilibrium population size, sigma^2 the
mutational variance, and s the invasion fitness gradient. This framework handles
frequency-dependent selection naturally and could model the coevolution of
reflector geometry (flower trait) and call frequency (bat trait), provided a
fitness function linking echo conspicuousness to pollination success can be
specified.

- Dieckmann, U. and Law, R. (1996). The dynamical theory of coevolution: a
  derivation from stochastic ecological processes. *Journal of Mathematical
  Biology*, 34, 579-612.
  DOI: [10.1007/BF02409751](https://doi.org/10.1007/BF02409751)

### 4.2 Mutualistic network trait evolution

Guimaraes et al. (2017) modelled trait evolution in mutualistic networks,
showing that network structure (nestedness, modularity) shapes the tempo and
mode of trait evolution. Their framework uses coupled Ornstein-Uhlenbeck
processes on a bipartite interaction network. Applied to bat-flower systems,
this would predict whether acoustic traits evolve towards a single optimum
(convergence) or diversify into distinct acoustic niches (divergence).

- Guimaraes, P. R., Pires, M. M., Jordano, P., Bascompte, J. and Thompson,
  J. N. (2017). Indirect effects drive coevolution in mutualistic networks.
  *Nature*, 550, 511-514.
  DOI: [10.1038/nature24273](https://doi.org/10.1038/nature24273)

### 4.3 Quantitative genetics with Gaussian matching

Nuismer, Jordano and Bascompte (2013) developed a quantitative genetics model
of coevolution in mutualistic networks where interaction strength depends on
trait matching via a Gaussian kernel:

> alpha_ij = exp(-(z_i - z_j)^2 / (2 * sigma_alpha^2))

where z_i and z_j are the traits of species i and j, and sigma_alpha controls
the specificity of matching. For bat-flower acoustics, the traits would be call
frequency and reflector resonance frequency, and alpha_ij would represent
pollination effectiveness.

- Nuismer, S. L., Jordano, P. and Bascompte, J. (2013). Coevolution and the
  architecture of mutualistic networks. *Evolution*, 67, 338-354.
  DOI: [10.1111/j.1558-5646.2012.01801.x](https://doi.org/10.1111/j.1558-5646.2012.01801.x)

### 4.4 Honest signalling and Spence game

Sun, Leshowitz and Rychtar (2018) modelled plant-pollinator communication as a
signalling game in which flowers invest in costly signals (nectar, scent,
colour) to attract pollinators. Acoustic reflectors represent an alternative
signalling modality. The framework predicts conditions under which honest
signalling (reflector conspicuousness correlates with nectar reward) is
evolutionarily stable.

- Sun, J., Leshowitz, M. I. and Rychtar, J. (2018). The signalling game
  between plants and pollinators. *Scientific Reports*, 8, 6686.
  DOI: [10.1038/s41598-018-24779-0](https://doi.org/10.1038/s41598-018-24779-0)

### 4.5 Optimal foraging and marginal value theorem

Howell and Hartl (1980) applied optimal foraging theory to nectar-feeding bats,
predicting patch residence times from the marginal value theorem. Their model
predicts the rate of diminishing returns at a flower patch but does not
incorporate detection time as a function of acoustic conspicuousness. Linking
reflector target strength to detection probability (and hence to the travel-time
component of the foraging cycle) would close this gap.

- Howell, D. J. and Hartl, D. L. (1980). Optimal foraging in glossophagine
  bats: when to give up. *American Naturalist*, 115, 696-704.
  DOI: [10.1086/283592](https://doi.org/10.1086/283592)

### 4.6 Bipartite network analysis

Queiroz et al. (2021) and Sritongchuay et al. (2019) constructed bipartite
networks of bat-flower interactions using field observation data. These studies
quantify network-level properties (connectance, nestedness, specialisation) but
define interaction strength purely from visitation frequency. None incorporate
acoustic interaction strength (e.g., target strength, spectral overlap between
call and reflector response) as an edge weight. Doing so would allow the
network structure to be predicted, rather than merely described, from physical
first principles.

---

## 5 Biomimetic engineering infrastructure

Several computational tools developed for biomimetic sonar provide
infrastructure directly applicable to bat-flower acoustic modelling.

### 5.1 BatSLAM

Steckel and Peremans (2013) combined the RatSLAM navigation algorithm with a
biomimetic sonar front-end to create BatSLAM, a simultaneous localisation and
mapping system inspired by bat echolocation. The sonar front-end uses a
gammatone filterbank and introduces the concept of an echolocation room transfer
function (ERTF), analogous to head-related transfer functions in human hearing.

- Steckel, J. and Peremans, H. (2013). BatSLAM: simultaneous localization and
  mapping using biomimetic sonar. *PLoS ONE*, 8, e54076.
  DOI: [10.1371/journal.pone.0054076](https://doi.org/10.1371/journal.pone.0054076)

### 5.2 AHRTF formalism

Schillebeeckx et al. (2011) formalised the acoustic head-related transfer
function (AHRTF) for bats, providing a complete description of the directional
filtering imposed by the bat's pinnae and noseleaf on both emission and
reception. This formalism is essential for any model that seeks to predict the
echo signal at the bat's ears (rather than at an idealised point receiver) from
a floral reflector.

- Schillebeeckx, F., De Mey, F. and Peremans, H. (2011). Bio-inspired sonar
  antennae: understanding bat-ear geometry for acoustic imaging. *International
  Journal of Robotics Research*, 30, 975-987.
  DOI: [10.1177/0278364910380474](https://doi.org/10.1177/0278364910380474)

### 5.3 SonoTraceLab

Jansen and Steckel (2024) developed SonoTraceLab, an open-source ray-tracing
simulator for ultrasonic environments with Monte Carlo diffraction modelling.
The simulator is implemented in MATLAB/CUDA and has been validated against the
FEM results of Simon et al. (2020) for the *Marcgravia evenia* dish-leaf. This
makes it the most directly relevant computational tool for simulating
bat-flower acoustic interactions at scale.

- Jansen, M. and Steckel, J. (2024). SonoTraceLab: a ray tracing toolbox for
  ultrasonic environments. *IEEE Access*, 12, 69371-69382.
  DOI: [10.1109/ACCESS.2024.3401226](https://doi.org/10.1109/ACCESS.2024.3401226)

### 5.4 PFM transform

Zhang, Du and Chen (2024) introduced the parametric frequency-modulated (PFM)
transform, a time-frequency analysis tool designed for chirp-like signals. The
PFM transform provides higher time-frequency resolution than the short-time
Fourier transform for FM sweeps, making it suitable for analysing both bat calls
and their echoes.

- Zhang, Y., Du, S. and Chen, X. (2024). Parametric frequency-modulated
  transform for bat echolocation signal analysis. *Journal of the Acoustical
  Society of America*, 156, 2596-2605.
  DOI: [10.1121/10.0032394](https://doi.org/10.1121/10.0032394)

### 5.5 TV-AR vocal tract model

Zhong et al. (2025) modelled bat vocalisations using a time-varying
autoregressive (TV-AR) model of the vocal tract, providing a parametric
description of call generation that can be coupled to propagation and
scattering models. This completes the forward model chain: call generation
(TV-AR) followed by propagation, then reflector scattering, then cochlear processing
(SCAT).

- Zhong, M., Chen, Y., Wang, H. and Li, S. (2025). Time-varying
  autoregressive modelling of bat echolocation calls. *Frontiers in Zoology*,
  22, 17.
  DOI: [10.1186/s12983-025-00573-3](https://doi.org/10.1186/s12983-025-00573-3)

---

## 6 Identified gaps and opportunities

The following seven open problems are ordered from the most concrete
(computational acoustics) to the most integrative (coevolutionary modelling).
Each is tractable with existing data and mathematical tools.

### 6.1 BEM scattering computation for floral reflectors

**Current state.** Only one numerical model exists: Simon et al. (2020) applied
FEM to an idealised spherical cap. Real floral geometries (vexilla, bell
corollas, dish-leaves) have not been modelled.

**Available data.** 3D surface scans and moulded replicas of floral reflectors
exist in several laboratories. Micro-CT data can provide sub-millimetre
resolution meshes.

**Mathematical tools.** Boundary element methods (BEM) are the standard tool for
exterior acoustic scattering. Open-source BEM libraries (e.g., Bempp, OpenBEM)
handle the Helmholtz equation on arbitrary surface meshes. The He = 3-30
regime is well within BEM capabilities without requiring fast multipole
acceleration.

**Contribution.** First-principles prediction of floral echo signatures from
geometry alone, enabling comparison with measured echoes and providing the
forward model needed for shape optimisation (Gap 6.2).

**Target journal.** *Journal of the Acoustical Society of America* or
*Journal of Theoretical Biology*.

### 6.2 Shape optimisation against biologically meaningful objectives

**Current state.** No study has optimised floral reflector geometry against any
objective function. The only shape study is Simon et al. (2020), which
parameterised spherical caps but did not perform optimisation.

**Available data.** The BEM forward model from Gap 6.1 provides the objective
function evaluation. Biologically relevant objectives include: maximising target
strength integrated over a specified angular range; maximising echo-to-clutter
ratio; maximising the spectral distinctiveness of the echo (quantified by
Kullback-Leibler divergence from a background clutter distribution).

**Mathematical tools.** Adjoint methods for shape sensitivity analysis in BEM;
gradient-based optimisation (L-BFGS) on a parameterised surface; topology
optimisation via level-set methods. Bayesian optimisation is an alternative when
gradients are unavailable.

**Contribution.** Quantitative prediction of whether observed floral shapes are
acoustically optimal, sub-optimal, or Pareto-optimal with respect to competing
objectives (e.g., acoustic conspicuousness vs. water retention).

**Target journal.** *Proceedings of the Royal Society B* or *Journal of the
Royal Society Interface*.

### 6.3 Porous media acoustic model for cephalium absorption

**Current state.** Simon et al. (2023) measured the absorption spectrum of the
*Espostoa frutescens* cephalium but fitted no physical model.

**Available data.** Absorption spectra at 20-100 kHz; SEM images providing fibre
diameter (~30 micrometre) and packing density; trichome length (~15 mm).

**Mathematical tools.** The Delany-Bazley empirical model relates absorption
coefficient to flow resistivity and frequency for fibrous materials. For a more
physically rigorous approach, Biot poroelastic theory separates frame and fluid
contributions. Johnson-Champoux-Allard (JCA) models offer an intermediate level
of complexity with five measurable parameters (porosity, tortuosity, flow
resistivity, viscous and thermal characteristic lengths).

**Contribution.** First physical model of acoustic absorption in a biological
trichome structure; prediction of absorption spectra from microstructural
parameters; identification of whether trichome geometry is tuned to pollinator
call frequency.

**Target journal.** *Journal of Experimental Biology* or *Journal of the
Acoustical Society of America*.

### 6.4 Parametric characterisation of glossophagine FM sweeps

**Current state.** FM sweep models (hyperbolic, exponential, logarithmic) have
been developed and fitted for insectivorous species (Section 2.3) but not for
glossophagine bats.

**Available data.** Call recordings from *Glossophaga soricina* (Simon et al.,
2006), *Leptonycteris yerbabuenae* (Gonzalez-Terrazas et al., 2016), and
additional species in publicly accessible bat call libraries.

**Mathematical tools.** Nonlinear least-squares fitting of instantaneous
frequency trajectories; information-theoretic model comparison (AIC, BIC) to
determine which functional form best describes glossophagine FM sweeps; Cramer-
Rao bound analysis to determine the theoretical range-resolution limits of each
waveform class.

**Contribution.** First systematic parameterisation of glossophagine FM sweeps;
determination of whether glossophagine sweeps are closer to Doppler-tolerant
(hyperbolic) or range-resolving (linear FM) designs; input parameters for the
SCAT model and for the POMDP observation function (Gap 6.6).

**Target journal.** *Journal of Comparative Physiology A* or *Journal of the
Acoustical Society of America*.

### 6.5 Adaptive dynamics model of acoustic trait coevolution

**Current state.** No model exists (Section 4). The component pieces are
available: BEM-derived fitness landscapes from reflector geometry (Gap 6.1),
FM sweep parameterisation (Gap 6.4), and the adaptive dynamics framework
(Section 4.1).

**Available data.** Phylogenetic trees for glossophagine bats and their host
plants; call frequency data across species; reflector morphometry. These
constrain initial conditions and allow model predictions to be tested against
observed trait distributions.

**Mathematical tools.** The canonical equation of adaptive dynamics (Dieckmann
and Law, 1996) for two coevolving traits (call centre frequency, reflector
curvature). Fitness is defined as pollination success, which depends on echo
signal-to-noise ratio, which in turn depends on the acoustic channel (call
spectrum convolved with reflector transfer function). Bifurcation analysis
identifies conditions for evolutionary branching (niche partitioning) vs.
convergence (trait matching).

**Contribution.** First mathematical model of bat-flower acoustic coevolution;
predictions testable against comparative phylogenetic data; identification of
whether acoustic traits are in evolutionary equilibrium or undergoing directional
selection.

**Target journal.** *American Naturalist* or *Evolution*.

### 6.6 POMDP formulation of flower approach

**Current state.** No POMDP model exists for glossophagine flower approach
(Section 2.6). The closest related work is the proportional navigation model
for insect interception (Ghose and Moss, 2006) and the neural network docking
model of Nguyen and Vanderelst (2025).

**Available data.** Flight trajectory data from high-speed video of flower
visits; call emission timing from on-board microphones; reflector echo
measurements (Sections 1.1-1.4).

**Mathematical tools.** The state space consists of 3D position and velocity of
the bat relative to the flower. The observation model maps the state to a
probability distribution over echo features (delay, spectrum, amplitude)
conditioned on the call parameters and reflector transfer function. The action
space includes flight control (3D acceleration) and call parameters (emission
timing, frequency range). Solving the POMDP yields an optimal policy that
specifies both flight manoeuvres and call adaptations as a function of the
current belief state. Point-based solvers (SARSOP, PBVI) handle the continuous
state space.

**Contribution.** Principled framework unifying perception and control in
flower approach; prediction of optimal call adaptation strategies; testable
predictions about whether bats achieve near-optimal performance.

**Target journal.** *PLoS Computational Biology* or *Journal of the Royal
Society Interface*.

### 6.7 Scattering-transform analysis of flower echoes

**Current state.** Echo classification has used SVMs on spectrograms (Yovel
et al., 2008) and CNNs (Simon et al., 2021) but not scattering transforms
(Section 3.4).

**Available data.** Echo recordings from the studies cited in Sections 1 and 3;
synthetic echoes generated by SonoTraceLab (Section 5.3).

**Mathematical tools.** The scattering transform (Mallat, 2012) applies iterated
wavelet decomposition followed by a modulus nonlinearity to produce
representations that are locally translation-invariant and stable to
deformations. For 1D signals (echo waveforms), the first-order scattering
coefficients correspond to a mel-like filterbank; the second-order coefficients
capture amplitude modulation structure. The log-frequency shift invariance
observed in bat echo perception (Section 3.2) maps naturally onto the
translation invariance of the scattering transform along the log-frequency axis.

**Contribution.** A biologically motivated, mathematically principled feature
extraction method for flower echoes; potentially superior classification
performance due to built-in invariances; formal connection between the
scattering transform's mathematical properties and the psychophysical invariance
observed in bat echo perception.

**Target journal.** *PLoS Computational Biology* or *IEEE Transactions on Signal
Processing*.

---

## References (alphabetical)

1. Altes, R. A. and Titlebaum, E. L. (1970). Bat signals as optimally Doppler
   tolerant waveforms. *Journal of the Acoustical Society of America*, 48,
   1014-1020. DOI: [10.1121/1.1912222](https://doi.org/10.1121/1.1912222)

2. Boonman, A. and Schnitzler, H.-U. (2005). Frequency modulation patterns in
   the echolocation signals of two vespertilionid bats. *Journal of Comparative
   Physiology A*, 191, 13-21.
   DOI: [10.1007/s00359-004-0566-8](https://doi.org/10.1007/s00359-004-0566-8)

3. Dieckmann, U. and Law, R. (1996). The dynamical theory of coevolution: a
   derivation from stochastic ecological processes. *Journal of Mathematical
   Biology*, 34, 579-612.
   DOI: [10.1007/BF02409751](https://doi.org/10.1007/BF02409751)

4. Ghose, K. and Moss, C. F. (2006). Steering by hearing: a bat's acoustic
   gaze is linked to its flight motor output by a delayed, adaptive linear
   law. *Journal of Neuroscience*, 26, 1704-1710.
   DOI: [10.1523/JNEUROSCI.4315-05.2006](https://doi.org/10.1523/JNEUROSCI.4315-05.2006)

5. Ghose, K., Horiuchi, T. K., Krishnaprasad, P. S. and Moss, C. F. (2006).
   Echolocating bats use a nearly time-optimal strategy to intercept prey.
   *PLoS Biology*, 4, e108.
   DOI: [10.1371/journal.pbio.0040108](https://doi.org/10.1371/journal.pbio.0040108)

6. Gonzalez-Terrazas, T. P., Koblitz, J. C., Fleming, T. H., Medellein, R. A.,
   Kalko, E. K. V., Schnitzler, H.-U. and Tschapka, M. (2016). How
   nectar-feeding bats localize their food. *PLoS ONE*, 11, e0163492.
   DOI: [10.1371/journal.pone.0163492](https://doi.org/10.1371/journal.pone.0163492)

7. Guimaraes, P. R., Pires, M. M., Jordano, P., Bascompte, J. and Thompson,
   J. N. (2017). Indirect effects drive coevolution in mutualistic networks.
   *Nature*, 550, 511-514.
   DOI: [10.1038/nature24273](https://doi.org/10.1038/nature24273)

8. Howell, D. J. and Hartl, D. L. (1980). Optimal foraging in glossophagine
   bats: when to give up. *American Naturalist*, 115, 696-704.
   DOI: [10.1086/283592](https://doi.org/10.1086/283592)

9. Jansen, M. and Steckel, J. (2024). SonoTraceLab: a ray tracing toolbox for
   ultrasonic environments. *IEEE Access*, 12, 69371-69382.
   DOI: [10.1109/ACCESS.2024.3401226](https://doi.org/10.1109/ACCESS.2024.3401226)

10. Masters, W. M., Jacobs, S. C. and Simmons, J. A. (1991). The structure of
    echolocation sounds used by the big brown bat *Eptesicus fuscus*. *Journal
    of the Acoustical Society of America*, 89, 1402-1413.

11. Ming, C., Haro, S., Bhatt, A. and Bhatt, D. (2021). A neural network
    model for echo processing in the inferior colliculus. *PLoS Computational
    Biology*, 17, e1008677.
    DOI: [10.1371/journal.pcbi.1008677](https://doi.org/10.1371/journal.pcbi.1008677)

12. Nuismer, S. L., Jordano, P. and Bascompte, J. (2013). Coevolution and the
    architecture of mutualistic networks. *Evolution*, 67, 338-354.
    DOI: [10.1111/j.1558-5646.2012.01801.x](https://doi.org/10.1111/j.1558-5646.2012.01801.x)

13. Saillant, P. A., Simmons, J. A., Dear, S. P. and McMullen, T. A. (1993).
    A computational model of echo processing. *Journal of the Acoustical
    Society of America*, 94, 2691-2712.
    DOI: [10.1121/1.407353](https://doi.org/10.1121/1.407353)

14. Schillebeeckx, F., De Mey, F. and Peremans, H. (2011). Bio-inspired sonar
    antennae. *International Journal of Robotics Research*, 30, 975-987.
    DOI: [10.1177/0278364910380474](https://doi.org/10.1177/0278364910380474)

15. Simon, R., Holderied, M. W., Koch, C. U. and von Helversen, O. (2011).
    Floral acoustics: conspicuous echoes of a dish-shaped leaf attract bat
    pollinators. *Science*, 333, 631-633.
    DOI: [10.1126/science.1204210](https://doi.org/10.1126/science.1204210)

16. Simon, R., Holderied, M. W. and von Helversen, O. (2006). Size
    discrimination of hollow hemispheres by echolocation in a nectar-feeding
    bat. *Journal of Experimental Biology*, 209, 3599-3609.
    DOI: [10.1242/jeb.02398](https://doi.org/10.1242/jeb.02398)

17. Simon, R., Rupitsch, S. J. and Baumann, M. J. (2020). Numerical simulation
    of bat-pollinated flower reflectors. *Proceedings of the National Academy of
    Sciences*, 117, 1367-1374.
    DOI: [10.1073/pnas.1909890117](https://doi.org/10.1073/pnas.1909890117)

18. Simon, R., Rupitsch, S. J., Guenther, F. and Kuc, R. (2021). Echo-based
    plant classification by echolocating bats. *PLoS Computational Biology*,
    17, e1009706.
    DOI: [10.1371/journal.pcbi.1009706](https://doi.org/10.1371/journal.pcbi.1009706)

19. Simon, R., Knornschild, M., Tschapka, M., Grund, A., Kalko, E. K. V. and
    Hoelzle, L. E. (2023). Acoustic absorption by cactus trichomes facilitates
    bat pollination. *Journal of Experimental Biology*, 226, jeb245263.
    DOI: [10.1242/jeb.245263](https://doi.org/10.1242/jeb.245263)

20. Steckel, J. and Peremans, H. (2013). BatSLAM: simultaneous localization
    and mapping using biomimetic sonar. *PLoS ONE*, 8, e54076.
    DOI: [10.1371/journal.pone.0054076](https://doi.org/10.1371/journal.pone.0054076)

21. Sun, J., Leshowitz, M. I. and Rychtar, J. (2018). The signalling game
    between plants and pollinators. *Scientific Reports*, 8, 6686.
    DOI: [10.1038/s41598-018-24779-0](https://doi.org/10.1038/s41598-018-24779-0)

22. Von Helversen, D. (2004). Object classification by echolocation in nectar
    feeding bats. *Journal of Comparative Physiology A*, 190, 515-521.
    DOI: [10.1007/s00359-004-0492-9](https://doi.org/10.1007/s00359-004-0492-9)

23. von Helversen, D. and von Helversen, O. (1999). Acoustic guide in
    bat-pollinated flower. *Nature*, 398, 759-760.
    DOI: [10.1038/19648](https://doi.org/10.1038/19648)

24. von Helversen, D. and von Helversen, O. (2003). Object recognition by
    echolocation: a nectar-feeding bat exploiting the flowers of a rain forest
    vine. *Journal of Comparative Physiology A*, 189, 327-336.
    DOI: [10.1007/s00359-003-0405-3](https://doi.org/10.1007/s00359-003-0405-3)

25. von Helversen, D., Holderied, M. W. and von Helversen, O. (2003). Echoes
    of bat-pollinated bell-shaped flowers. *Journal of Experimental Biology*,
    206, 1025-1034.
    DOI: [10.1242/jeb.00203](https://doi.org/10.1242/jeb.00203)

26. Yovel, Y., Franz, M. O., Stilz, P. and Schnitzler, H.-U. (2008). Plant
    classification from bat-like echolocation signals. *PLoS Computational
    Biology*, 4, e1000032.
    DOI: [10.1371/journal.pcbi.1000032](https://doi.org/10.1371/journal.pcbi.1000032)

27. Zhang, Y., Du, S. and Chen, X. (2024). Parametric frequency-modulated
    transform for bat echolocation signal analysis. *Journal of the Acoustical
    Society of America*, 156, 2596-2605.
    DOI: [10.1121/10.0032394](https://doi.org/10.1121/10.0032394)

28. Zhong, M., Chen, Y., Wang, H. and Li, S. (2025). Time-varying
    autoregressive modelling of bat echolocation calls. *Frontiers in Zoology*,
    22, 17.
    DOI: [10.1186/s12983-025-00573-3](https://doi.org/10.1186/s12983-025-00573-3)
