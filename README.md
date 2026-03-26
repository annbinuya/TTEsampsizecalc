# Time-to-Event Sample Size Calculator

An interactive **R Shiny** application for designing **time-to-event studies** under **exponential survival assumptions**, using **log-rank** and **Schoenfeld / Lachin-Foulkes**-based formulas.

## Author

| Field | Details |
|---|---|
| Author | **Mary Ann Binuya** |

## Overview

This app supports rapid exploration of time-to-event trial designs through an interactive interface that links assumptions, formulas, and visual outputs in one place.

It is designed to help users move quickly from clinical assumptions to operational quantities such as required events, sample size, event probability, and achievable power under budget constraints.

## What the app can do

| Capability | Description |
|---|---|
| **Single-arm design calculations** | Calculates required events and sample size for one-arm time-to-event studies against a benchmark or historical control assumption. |
| **Two-arm design calculations** | Supports randomized two-arm designs with treatment and control inputs. |
| **Multiple hypothesis settings** | Includes **equality**, **superiority**, **noninferiority**, and **equivalence** options. |
| **Flexible input formats** | Accepts assumptions as **median survival**, **hazard rate**, or **landmark survival probability** at a user-defined time point. |
| **Hazard ratio derivation** | Converts inputs into a consistent hazard-based framework and displays the implied **HR** and log(HR). |
| **Required events calculation** | Estimates the number of events needed for the target alpha and power. |
| **Sample size estimation** | Converts event requirements into total sample size using accrual, follow-up, and dropout assumptions. |
| **Allocation ratio support** | For two-arm trials, allows unequal randomization through user-defined treatment:control allocation ratio. |
| **Dropout-adjusted design** | Incorporates anticipated dropout/loss to follow-up into effective event probability and sample size. |
| **Budget back-calculation** | Given a maximum feasible sample size, estimates the expected number of events, achieved power, and alpha relaxation needed to hit target power. |
| **Sensitivity analysis** | Evaluates how required events and sample size change across nearby hazard ratio scenarios around the selected design. |
| **Visual design diagnostics** | Provides plots for survival curves, test-statistic distributions, cumulative events, and a power heatmap. |
| **Transparent calculation steps** | Shows step-by-step derivations to make the design logic easier to review, explain, and present. |

## Key features

### 1. Flexible design input
Users can specify assumptions in whichever format is most natural for the setting:

| Input mode | Example use |
|---|---|
| **Median survival** | When protocol assumptions are expressed in median PFS or OS |
| **Hazard rate** | When a rate-based assumption is already available |
| **Landmark survival S(t\*)** | When survival probability at a specific time point is the main planning quantity |

### 2. Core design outputs
The app provides both statistical and operational outputs.

| Output | Description |
|---|---|
| **Required events (D)** | Event count needed to achieve the target design performance |
| **Sample size (N)** | Total required enrollment after accounting for event probability and dropout |
| **Arm-specific sample size** | For two-arm designs, total N is split into treatment and control counts |
| **Hazard ratio (HR)** | Derived from the user’s design assumptions |
| **Dropout-adjusted event probability** | Effective event probability used to translate D into N |
| **Achieved power** | Recalculated power after event rounding |
| **Budget-constrained operating characteristics** | Expected events, achievable power, and required alpha relaxation at a fixed N |

### 3. Sensitivity analysis
The app can explore how the design behaves when assumptions change modestly around the selected hazard ratio.

This helps users assess:

- how fragile or robust the planned design is
- whether small changes in treatment effect materially alter required sample size
- whether the proposed design remains practical across plausible scenarios

### 4. Budget-aware planning
The app includes a **budget back-calculation module**, which is useful when enrollment is limited by feasibility rather than by a purely statistical target.

This allows users to ask:

- What power can I achieve with my maximum feasible N?
- How many events would I expect under this constraint?
- How much would alpha need to change to recover the target power?

This is especially useful during early protocol discussions, internal planning, or feasibility negotiations.

### 5. Visual outputs
The app includes several visual summaries to make assumptions and consequences easier to understand.

| Visualization | Purpose |
|---|---|
| **Survival curves** | Shows the implied exponential survival functions under the chosen assumptions |
| **Test-statistic distribution** | Illustrates the null and alternative distributions and the relationship between alpha, beta, and power |
| **Cumulative events plot** | Shows event accumulation over calendar time under accrual and follow-up assumptions |
| **Power heatmap** | Displays how power changes across sample size and hazard ratio combinations |

These visuals can help with both technical review and communication with collaborators who prefer a more intuitive representation of the design.

## Statistical basis

The app is based on standard asymptotic methods for time-to-event design under exponential survival assumptions.

| Methodological component | Basis |
|---|---|
| Required events for log-rank type designs | **Schoenfeld (1981, 1983)** |
| Event probability under accrual and follow-up | **Lachin and Foulkes (1986)** |
| Exponential survival parameterization | Conversion between median, hazard rate, and survival probability |

## Version history

| Version | Date | Summary |
|---|---|---|
| **v1** | **2025-10-08** | Initial **one-arm** version |
| **v2** | **2026-03-26** | Added **two-arm** functionality; formulas for some two-arm settings remain work in progress (DO NOT USE!) |

## Intended use

This app is intended for:

| Use case | Description |
|---|---|
| **Early trial design exploration** | Rapidly compare assumptions and operating characteristics |
| **Sensitivity assessment** | Evaluate how design requirements change under nearby effect sizes |
| **Feasibility discussions** | Examine trade-offs between statistical targets and practical enrollment limits |
| **Educational use** | Demonstrate how survival assumptions map to events, sample size, and power |
| **Protocol support** | Help structure and communicate design assumptions before formal statistical finalization |

It is particularly relevant for oncology and other clinical trial settings involving endpoints such as:

- **PFS**
- **OS**
- **DFS**
- **EFS**
- other time-to-event endpoints

## Requirements

| Package | Purpose |
|---|---|
| `shiny` | Main application framework |
| `shinyjs` | UI interactivity and conditional display logic |

Install with:

```r
install.packages(c("shiny", "shinyjs"))
