---
title: "Survival in Seattle: Earthquake recurrence intervals and time-dependent 
earthquake likelihoods from non-parametric survival analysis of 
paleoearthquakes in the Puget Lowland (WA, USA) and San Andreas Fault (CA, 
USA)"

author:
- name: Richard Styron
  affiliation: Earth Analysis
  email: richard.h.styron@gmail.com
  
- name: Kate Scharer
  affiliation: US Geological Survey
  
- name: Brian Sherrod
  affiliation: US Geological Survey
  
abstract: "abstract"
...


# Introduction

The characterization of earthquake occurrence in time is of obvious importance 
in seismic hazard analysis, and is also of substantial (though less immediately 
practical) interest in the insight that it yields into the physics of the 
earthquake process: Different types of earthquake recurrence behaviour imply 
different types of loading, triggering and unloading (**reword**) of faults or 
fault systems [e.g. @faenza_non-parametric_2003] (**say something about 
interaction perhaps**). These 


## Recurrence models

![Map of fault traces in the Puget lowlands hosting the earthquakes analyzed 
here. \label{fig_eq_map}](./figures/rupture_map.pdf)


# Earthquake geology of the Puget Lowlands and San Andreas

## Puget Lowlands

## San Andreas

# Empirical recurrence model construction

We construct empirical earthquake recurrence models in a straightforward 
manner:

1. Draw a set of $N$ random samples directly from each of the paleoearthquake 
   age PDFs through an inverse transform sampling algorithm.

2. Make $N$ sequences of earthquakes by ordering one sample per earthquake from 
   each of the earthquake age sample sets (the nature of this ordering depends 
   on availability of relative time constraints and will be discussed below).

3. Calculate the interevent times (i.e. the time interval between successive 
   events) for all of the earthquakes in each of the $N$ sequences

4. Concatenate the interevent times, sort them, and create a recurrence 
   (interevent) PDF through a kernel density estimation of the interevent 
   times.


### Sorting vs culling based on geologic relationships

An important consideration here is the presence of stratigraphic constraints on 
paleoearthquake ordering when earthquake age PDFs overlap. If the earthquakes 
have no stratigraphic or other geologic ordering constraints (e.g., they are 
from different trenches), each sequence of earthquakes can simply be sorted and 
then differenced to obtain the interevent time sequence. However, if the ages 
overlap but Earthquake A is demonstrably older than the next (for example 
through superposition of syn-earthquake colluvial wedges), then the sample age 
sequence cannot simply be sorted, because some sample earthquake sequences 
would place the date for the older Earthquake A *after* the age for the younger 
Earthquake B. In this situation, we instead reject sample earthquake sequences 
that are not in geologic order.

In the Puget Lowland, this is not a concern: The earthquakes from each 
individual fault trace are well-separated in time, and no stratigraphic or 
cross-cutting relationships are known to constrain relative ages of temporally 
overlapping events on different faults, so these can be treated independently.

However, geologic (mostly stratigraphic) relative age constraints are strong in 
the San Andreas data. For these analyses, we cull the sample sets that violate 
relative age constraints.

## Puget Lowlands recurrence

\begin{figure}[tb]
  \centering
  \includegraphics{./figures/pug_recurrence.pdf}
  \caption{Earthquake timing and recurrence intervals in the Puget Lowlands. 
  \textbf{a}: Time probabilities for each earthquake from the OxCal modeling. 
  \textbf{b}: Interevent times for each set of consecutive earthquakes (thin 
  lines) and kernel-density estimate for all interevent times (thick line). 
  Note that the left \emph{y}-axis scale corresponds to the individual 
  interevent PDFs, while the right \emph{y}-axis scale corresponds to the total 
  kernel density estimated recurrence PDF.\label{fig_pug_recurrence}}
  \end{figure}

Earthquakes in the Puget Lowland are dispersed throughout the past 16 ka, with 
most of the earthquakes happening in the past 4000 years (Figure 
\ref{fig_pug_recurrence}a). Given the considerable overlap in the age PDFs 
(especially for the older events), the earthquakes appear to be grouped in 
time, with the oldest group at 14-10 ka, another group at 9-6 ka, and a final 
group from 4 ka to the present. 


### Seattle Fault Zone

\begin{figure}[tb]
  \centering
  \includegraphics{./figures/sfz_recurrence.pdf}
  \caption{Earthquake timing and recurrence intervals in the Seattle Fault 
  Zone. See Figure \ref{fig_pug_recurrence} for explanation. 
  \label{fig_sfz_recurrence}}
\end{figure}


The Seattle Fault Zone has a billion earthquakes per year (Figure 
\ref{fig_sfz_recurrence}.

##  San Andreas recurrence

# Survival analysis

Survival analysis is a branch of statistics that deals with the time to events, 
commonly deaths (hence its name) or product failures. It is very commonly used 
in the life sciences (especially medicine and public health) as well as in 
engineering (where it is often called 'reliability analysis'), and has seen 
application in a wide variety of fields. Unfortunately, this has lead to some 
variability in terminology that can make searching for certain ideas 
challenging; we will provide a list of common terms for the same concept where 
appropriate.

The techniques of survival analysis are well-suited for studying a range of 
geoscientific phenomena, but there seems to be little explicit mention of it in 
the literature or the classroom. Several of the key functions (particularly the 
hazard function) are used commonly enough in the earthquake science community 
[e.g., @davis_longer_1989; @reasenberg_earthquake_1989; @sornette_paradox_1997; 
@matthews_brownian_2002] that the terminology is not universally unfamiliar, 
but in our science, these tools seem to be used almost exclusively 
by a small set of specialists in statistical seismology who are already well 
versed in them. As a consequence, it is hard to find a straightforward 
explanation of the principles and general description of the key relationships 
in the geoscience literature.

Below, we give a simple overview of the aspects of survival analysis used in 
this paper. This should be considered an extremely superficial introduction to 
a mature science, but (as is often the case with applied mathematics) the 
basics can be quite helpful in both clarifying our ideas and giving us powerful 
quantitative tools. Additional techniques in survival analysis can be used if 
an application calls for it; for example -@faenza_non-parametric_2003 apply a 
proportional hazard model to a study of Italian seismicity to incorporate the 
effects of earthquake magnitude and spatial occurrence into a time- and 
space-dependent model.

As we are operating on empirical PDFs that have no simple analytical form, the 
mathematical descriptions will be general, which should aid in comprehension, 
perhaps at the expense of mathematical elegance.

## Survival analysis mathematics

Survival analysis is based on a few simple equations incorporating the 
recurrence interval PDF, or $f(t)$, which is the probability that some event 
will occur at time $t$ from a start time (such as a birth; for our purposes, 
this is the time since the last event). First, we define $T \ge 0$, a 
non-negative random variable representing the time between earthquakes.

Next, we can define $F(t)$, which is the probability that $T$ is as short or 
shorter than some $t$, i.e. that a child survives for fewer than $t$ years or 
that the time between earthquakes is less than $t$ years:
\begin{equation}
  \label{eqn_F_def}
  F(t) = \begin{cases} \Pr[T \le t] & t \ge 0 \\ 0 & t < 0 \end{cases}
\end{equation}
$F(t)$ is the *cumulative distribution function* of $f(t)$, and is easily 
calculated from it if $f(t)$ is known:
\begin{equation}
  \label{eqn_F_int}
  F(t) = \int_{-\infty}^t f(t) dt \;.
\end{equation}
In the case of empirical (or non-analytical) $f(t)$, as in this study, $F(t)$ 
can be calculated through numerical integration.

Then, we can define the *survival function* (also known as the reliability 
function) $S(t)$ as
\begin{equation}
  \label{eqn_survival_fn}
  S(t) = 1 - F(t) \; .
\end{equation}
The survival function, as the compliment of $F(t)$, describes the probability 
that $T$ is longer than $t$, or $\Pr[T>t]$, i.e. that a child will live beyond 
a certain age. If $f(t)$ or $F(t)$ is not known but some some samples of $T$ 
are known, $S(t)$ can still be estimated, for example through the Kaplan-Meier 
Estimator [@kaplan_nonparametric_1958] (which can also handle *censored data*, 
discussed below). Though the survival function is foundational, we will not use 
it directly in this work. Instead, several additional functions are more 
relevant for seismic hazard.

A more immediately useful function, the *hazard function* (also called the 
*hazard rate function* or *failure function*) $\lambda (t)$ describes the 
instantaneous probability of an event at time $t$ given that the event has not 
yet occurred. $\lambda(t)$ can be derived from the previous functions as:
\begin{equation}
  \label{eqn_haz}
    \lambda(t) = f(t) / S(t) \;.
\end{equation}
Assuming that the last date of an earthquake on a fault is known, $\lambda(t)$ 
may be used in probabilistic earthquake hazard analysis. 

A plot of $t$ vs. $\lambda(t)$ may assume many forms, and these reflect the 
processes controlling the distribution of lifetimes or interevent times, or the 
assumptions thereof. 

### Censored data

### operations on KDE functions

## Puget hazard

### SFZ hazard

## SAF hazard

- wrightwoood

- pallet creek

# Discussion

# Conclusions

# References
