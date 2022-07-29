+++
title = "Minilab 2: Future evolution of HMXBs"
hascode = true
date = Date(2019, 3, 22)
rss = "Model"
+++
@def tags = ["syntax", "code"]

# Minilab 2: Future evolution of HMXBs

\toc

## High mass X-ray binaries with black holes

The objective of this lab will be to study the future evolution of X-ray binaries known 
to have black holes, for which we also have a well characterised orbital solution. The following
table shows the properties of these systems that will be relevant to our models:

| Name  | $M_\mathrm{BH}$ $[M_\odot]$ | $M_\mathrm{star}$ $[M_\odot]$ | Period $[d]$ | Reference |
|-------|-------------------------------|-----------------------------|----------------|---|
| Cyg X-1 | $21.2^{+2.2}_{-2.3}$ | $40.6^{+7.7}_{-7.1}$| $5.60$ | \cite{Miller-Jones+21}
| LMC X-1 | $10.9\pm1.4$ | $31.79\pm3.67$| $3.91$ | \cite{Orosz+09} |
| M33 X-7 | $15.65\pm1.45$ | $70.0\pm6.9$ | $3.45$ | \cite{Orosz+07}|
\tablecaption{HMXBs with orbital solutions that are known to host black holes.}

All of these systems are expected to have started as a non-degenerate binary, where the
progenitor of the black hole transferred mass to the star that is currently observed in the
system. This potentially means that the star that is observed digresses from the evolution
of a single star, and if we wanted to be fully consistent we would need to model the evolution
of the system starting from two stars at the zero age main sequence. Here we will avoid this
complication by approximating the star as a zero age main sequence star of the present-day mass, and we will also
ignore the errors in the observations. So, for example, to model Cyg X-1 we will use the following
in our `inlist_project`:

```fortran
m1 = 40.6d0
m2 = 21.2d0
initial_period_in_days = 5.6d0
```
The objective will be to model the three systems in the table above, study how binary interaction differs between them,
and assess some of the limitations of the `binary` module.

## Modifications to the template
We will start this minilab using the solution to minilab 1. Two small adjustements that we will need are the following
- Set `max_model_number=400` in `inlist1`. This is enough to model the systems until near-detachment. We will ignore the core helium burning phase.
- Add `Kipp_xaxis_name = 'star_age'` and change `History_Panels1_xaxis_name='star_age'` in `inlist1`. Having properties shown against time instead of model number in `pgstar` will help understand the different timescales that come into play.

@@colbox-blue 
Did you have issues finishing minilab1? You can get a working copy of its solution [here](/assets/minilab2/template2.tar.gz).
@@

## Run you models
Now, using the adjusted template, run all three simulations. Answer the following questions, and discuss them with those on your table. You can also
spread out the work among your nearby colleagues, if so be sure to have two people modelling M33 X-7:
- How much does the mass transfer phase last compared to the total lifetime of the system (consider here only the time modeled within the `max_model_number=400` constraint).
- Do the evolution of these systems differ qualitatively from the one on minilab1?
- Are all simulations running smoothly? Or are you running into numerical issues?

The last point in particular will be dealt with in the next section, so if you find everything is running smoothly be sure to double check!

## Limitations of MESA

As you might have seen, one of the systems modeled runs into some issues (to put it lightly). Mass transfer rates go to very high values, approaching
a solar mass per year. This points to some limitations within MESA, as the code models the donor star as a 1-dimensional object. For extreme mass ratios,
the shrinkage of the orbit as mass transfer proceeds becomes extreme enough that the star cannot adjust itself through mass loss to avoid extreme overflow.
In such a case the donor would very likely engulf its companion, initiating a process of common-envelope evolution which is fundamentally 3-dimensional.
MESA being a 1-dimensional code cannot deal with such a situation, but rather tries to keep modeling this evolutionary phase as a stable mass transfer event
with an ever increasing mass transfer rate, which eventually leads to numerical problems.

Although there are ways to approximate a common-envelope phase with MESA, here we wish to simply construct a physical way to identify when an unstable
mass transfer phase could start, and terminate the evolution at that stage.

## References

* \biblabel{Miller-Jones+21}{Miller-Jones et al. (2021)} [Miller-Jones, James C. A., et al.](https://ui.adsabs.harvard.edu/abs/2021Sci...371.1046M/abstract), Science, Volume 371, Issue 6533, pp. 1046-1049 (2021)
* \biblabel{Orosz+07}{Orosz et al. (2007)} [Orosz, Jerome A., et al.](https://ui.adsabs.harvard.edu/abs/2007Natur.449..872O/abstract), Nature, Volume 449, Issue 7164, pp. 872-875 (2007)
* \biblabel{Orosz+09}{Orosz et al. (2009)} [Orosz, Jerome A., et al.](https://ui.adsabs.harvard.edu/abs/2009ApJ...697..573O/abstract), The Astrophysical Journal, Volume 697, Issue 1, pp. 573-591 (2009)