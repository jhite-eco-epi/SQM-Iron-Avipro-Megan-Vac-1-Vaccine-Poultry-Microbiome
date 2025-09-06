# Beta Diversity: Jaccard (NMDS + PERMANOVA)

This folder contains Jaccard β-diversity analyses for the poultry vaccination × iron-supplement experiment.

* **Script:** `stats_jaccard.R`

## PERMANOVA (adonis2)

| Effect | df | Sum Sq | R² | *F* | *p*-value |
|--------|:--:|-------:|----:|----:|----------:|
| Vaccination | 1 | 0.3414 | 0.054 | 1.73 | **0.002** |
| Supplement  | 2 | 0.7985 | 0.126 | 2.03 | **0.001** |
| Vaccination × Supplement | 2 | 0.6749 | 0.106 | 1.71 | **0.001** |
| Residual | 23 | 4.5336 | 0.714 | — | — |
| **Total** | 28 | 6.3485 | 1.000 | — | — |

### Interpretation

* Both vaccination and iron supplementation independently explain significant variation in community presence/absence patterns.
* A significant interaction suggests the effect of supplementation differs between vaccinated and unvaccinated birds.

## Pair-wise PERMANOVA (Bonferroni-adjusted)

None of the 15 pair-wise comparisons retained significance after Bonferroni correction (*adj. p* > 0.05). The largest F-ratios involved contrasts between Control vs FeSO₄ within the same vaccination status, but they did not survive multiple-testing correction.

## NMDS Visualization

Non-metric multidimensional scaling (NMDS) ordination was used to visualize the Jaccard β-diversity patterns identified through PERMANOVA analysis.

### Dimensionality Optimization for Visualization

To determine the optimal number of dimensions for NMDS visualization, stress values were tested across 1–6 dimensions using `nmds_dimensionality_jaccard.R`:

| Dimensions | Stress | Converged | Stress Reduction | % Reduction |
|:----------:|-------:|:---------:|-----------------:|------------:|
| 1D | 0.2825 | Yes | — | — |
| 2D | 0.1474 | Yes | 0.1351 | 47.8% |
| 3D | 0.1150 | Yes | 0.0324 | 22.0% |
| **4D** | **0.0812** | **Yes** | **0.0338** | **29.4%** |
| 5D | 0.0661 | Yes | 0.0151 | 18.6% |
| 6D | 0.0532 | Yes | 0.0129 | 19.5% |

**Selected: 4D NMDS** (stress = 0.081) provides optimal representation quality with stress < 0.10 (good fit). For visualization purposes, the first two dimensions of the 4D solution are plotted with treatment groups distinguished by color and shape, complemented by 95% confidence ellipses.

### Visualization Scripts

* **Dimensionality analysis:** `nmds_dimensionality_jaccard.R`
* **NMDS plotting:** `nmds_plot_jaccard.R`
* **Composite figures:** `composite_nmds_stress_analysis.R`, `composite_beta_diversity_plot.R`

---
For detailed outputs see comments embedded in `stats_jaccard.R`.
