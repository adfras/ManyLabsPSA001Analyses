# Intro5: Research Aims, Questions, and Hypotheses (Location-Scale Variability)

## Purpose (Research Aims)
This project examines how trial-level variability and sampling across random facets (participants, stimuli, and sites) shape the interpretation of psychological effects. Using hierarchical location-scale models fit at the trial level, we estimate mean effects and within-person volatility jointly while allowing effects to vary across people and, where relevant, across labs and stimuli. The objective is diagnostic rather than corrective: to determine whether apparent instability reflects mean shifts, variance shifts, or differences in precision across facets.

## Operationalization and evaluation
We define within-person volatility as the participant-specific residual standard deviation estimated by the scale component of the hierarchical model. Directional stability is summarized using the posterior probability that an individual's slope is positive; individuals are classified as directionally positive when P(slope > 0) exceeds a pre-specified threshold (e.g., 0.95), and as sign-uncertain otherwise. To test whether modeled variability reflects genuine structure rather than overfitting, models are compared using held-out trials with a proper scoring rule (expected log predictive density), alongside calibration checks such as predictive interval coverage. The primary comparison is between a location-scale model and a homoskedastic mixed-effects baseline that shares the same mean structure but constrains residual variance to be constant across individuals.

## Research Questions

**RQ1 (Task comparison):**
How do within-person volatility and person-level effect heterogeneity differ between Stroop and PSA001 under the same hierarchical location-scale framework?

**RQ2 (Direction / prevalence):**
For each task, what is the prevalence of positive, near-zero, and negative person-level effects once trial noise and person-specific variance are modeled?

**RQ3 (Held-out prediction):**
Does modeling person-specific residual variance improve out-of-sample predictive performance (and calibration) relative to a homoskedastic mixed-effects baseline, particularly in the more variable task?

## Hypotheses

**H1 (RQ1/RQ2):** Stroop will show a high prevalence of positive individual effects (high posterior probability of a positive slope) and a tighter person-level slope distribution than PSA001.

**H2 (RQ1/RQ2):** PSA001 will show higher within-person volatility and/or greater heterogeneity in volatility than Stroop, yielding more individuals with uncertain direction (greater sign-uncertainty) even if the group mean effect is present.

**H3 (RQ3):** The location-scale model will yield better-calibrated predictive distributions than a homoskedastic baseline on held-out trials (e.g., higher expected log predictive density and/or better predictive interval coverage), with larger gains in PSA001.
