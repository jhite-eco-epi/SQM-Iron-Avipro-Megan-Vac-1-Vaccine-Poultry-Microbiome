# Beta Diversity: Bray–Curtis (NMDS + PERMANOVA)

This folder contains Bray–Curtis β-diversity analyses for the poultry vaccination × iron-supplement experiment.

* **Script:** `stats_bray_curtis.R`

## PERMANOVA (adonis2)

| Effect | df | Sum Sq | R² | *F* | *p*-value |
|--------|:--:|-------:|----:|----:|----------:|
| Vaccination | 1 | 0.2328 | 0.071 | 2.65 | **0.001** |
| Supplement  | 2 | 0.6032 | 0.184 | 3.43 | **0.001** |
| Vaccination × Supplement | 2 | 0.4143 | 0.127 | 2.36 | **0.001** |
| Residual | 23 | 2.0221 | 0.618 | — | — |
| **Total** | 28 | 3.2723 | 1.000 | — | — |

### Interpretation

* Both vaccination and iron supplementation independently explain significant variation in community structure.
* A significant interaction suggests the effect of supplementation differs between vaccinated and unvaccinated birds.

## Pair-wise PERMANOVA (Bonferroni-adjusted)

Only the most notable adjusted *p*-values (< 0.05) are shown.

| Comparison | *F* | adj. *p* | Significance |
|------------|----:|---------:|:------------:|
| Vaccinated.Control vs Vaccinated.SQMFe | 5.67 | 0.045 | • |
| Unvaccinated.Control vs Unvaccinated.FeSO₄ | 10.87 | 0.075 | — |
| Vaccinated.Control vs Unvaccinated.SQMFe | 2.60 | 0.30 | — |

The single significant (•) comparison indicates SQMFe alters the vaccinated birds' microbiota relative to their control counterparts after multiple-test correction. Other contrasts did not reach significance once corrected.

## NMDS Visualization

Non-metric multidimensional scaling (NMDS) ordination was used to visualize the Bray-Curtis β-diversity patterns identified through PERMANOVA analysis.

### Dimensionality Optimization for Visualization

To determine the optimal number of dimensions for NMDS visualization, stress values were tested across 1–6 dimensions using `nmds_dimensionality_bray.R`:

| Dimensions | Stress | Converged | Stress Reduction | % Reduction |
|:----------:|-------:|:---------:|-----------------:|------------:|
| 1D | 0.3441 | Yes | — | — |
| 2D | 0.2058 | Yes | 0.1383 | 40.2% |
| 3D | 0.1327 | Yes | 0.0731 | 35.5% |
| **4D** | **0.0946** | **Yes** | **0.0381** | **28.7%** |
| 5D | 0.0718 | Yes | 0.0228 | 24.1% |
| 6D | 0.0550 | Yes | 0.0168 | 23.4% |

**Selected: 4D NMDS** (stress = 0.095) provides optimal representation quality with stress < 0.10 (good fit). For visualization purposes, the first two dimensions of the 4D solution are plotted with treatment groups distinguished by color and shape, complemented by 95% confidence ellipses.

### Visualization Scripts

* **Dimensionality analysis:** `nmds_dimensionality_bray.R`
* **NMDS plotting:** `nmds_plot_bray.R`  
* **Composite figures:** `composite_nmds_stress_analysis.R`, `composite_beta_diversity_plot.R`

---
For detailed outputs see comments embedded in `stats_bray_curtis.R`.
