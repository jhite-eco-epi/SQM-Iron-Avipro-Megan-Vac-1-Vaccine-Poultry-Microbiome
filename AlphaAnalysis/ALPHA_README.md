# Alpha Diversity Analysis

This directory contains alpha diversity analyses for the poultry vaccination × iron supplementation experiment, examining within-sample diversity using observed richness (ASVs) and Shannon diversity index.

## Key Files

* **Main Script:** `composite_alpha_diversity_plot.R`
* **Statistical Tables:** `observed_richness_stats.csv`, `shannon_diversity_stats.csv`
* **Individual Plots:** `shannon_by_treatment.png`, `shannon_diversity_plot.png`

The composite figure shows both alpha diversity metrics side by side:
- **Panel A:** Observed richness (amplicon sequence variants, ASVs)  
- **Panel B:** Shannon diversity index

Both panels display treatment combinations with vaccination status on the x-axis and iron supplementation indicated by fill color (Control, FeSO₄, SQM® Iron).

## Statistical Analysis

Scheirer-Ray-Hare tests (non-parametric two-way ANOVA analogue) were applied to examine the individual and interactive effects of vaccination and supplementation on both alpha diversity metrics.

### Observed Richness (ASVs)

| Effect | df | Sum Sq | H statistic | *p*-value |
|--------|:--:|-------:|------------:|-----------:|
| Vaccination | 1 | 0.21 | 0.003 | 0.957 |
| Supplement  | 2 | 106.98 | 1.477 | 0.478 |
| Vaccination × Supplement | 2 | 245.42 | 3.388 | 0.184 |

### Shannon Diversity Index

| Effect | df | Sum Sq | H statistic | *p*-value |
|--------|:--:|-------:|------------:|-----------:|
| Vaccination | 1 | 88.86 | 1.226 | 0.268 |
| Supplement  | 2 | 19.51 | 0.269 | 0.874 |
| Vaccination × Supplement | 2 | 8.74 | 0.121 | 0.942 |

### Interpretation

* Neither vaccination nor iron supplementation produced statistically significant main effects on either alpha diversity metric (all *p* > 0.05).
* The interaction terms were also non-significant for both metrics, indicating no evidence that the effects of supplementation depend on vaccination status.
* Both observed richness and Shannon diversity show similar patterns of non-significant treatment effects.
