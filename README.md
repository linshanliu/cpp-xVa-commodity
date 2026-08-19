# Commodity XVA Engine (Monte Carlo)

A portfolio project simulating the core quantitative workflow of a commodities XVA / exotics desk:
pricing counterparty credit and funding risk (CVA/FVA) on a netting set of commodity derivatives,
using Monte Carlo exposure simulation under a two-factor stochastic convenience-yield model.

This project is a companion to two earlier projects:
- **FDM PDE Pricer** (5-component architecture: PDE / Payoff / BoundaryCondition / Grid / Scheme)
- **SVI Volatility Surface Calibrator** (arbitrage-checked, tested on real options data)

Where those projects covered low-dimensional PDE pricing and vol surface construction, this project
covers the high-dimensional simulation side of the desk: exposure profiles, netting, collateral, and
XVA aggregation — the reason Monte Carlo (not PDE) is the workhorse method for this kind of book.

---

## Motivation

Commodities XVA desks price the credit and funding cost of a *portfolio* of trades against a single
counterparty, not a single trade in isolation. Because the relevant quantity — the netting set's
future mark-to-market distribution across many trades, tenors, and (often correlated) underlyings —
is high-dimensional and path-dependent, PDE methods become impractical and Monte Carlo becomes the
natural tool. This project builds a simplified but structurally faithful version of that pipeline, to
build genuine intuition (not just textbook familiarity) for concepts covered in commodities XVA
interviews: exposure metrics, netting, collateral (CSA) mechanics, wrong-way risk, and the numerical
challenges of simulating and aggregating exposure profiles.

## What this project will actually do

1. **Simulate** commodity spot/convenience-yield dynamics under a risk-neutral measure using a
   two-factor Gibson-Schwartz model (mean-reverting convenience yield + spot), discretized via
   Euler or Milstein, with correlation across a small number of commodities (e.g. crude + nat gas)
   via Cholesky decomposition.
2. **Reprice** a small netting set (3-5 trades: forwards, European and American/swing-style options)
   along every simulated path at every time step, using Longstaff-Schwartz regression for
   early-exercise instruments.
3. **Aggregate** exposure into standard metrics: Expected Exposure (EE), Potential Future Exposure
   (PFE), Expected Positive/Negative Exposure (EPE/ENE), with and without a simplified CSA
   (threshold, minimum transfer amount, margin period of risk).
4. **Compute XVA**: CVA from the EE profile against a bootstrapped counterparty default curve; FVA
   from a funding spread over the exposure profile; a simplified wrong-way risk overlay linking
   counterparty hazard rate to the commodity price path.
5. **Address the numerical/engineering problems** a real desk faces at smaller scale: memory
   management for large path x time-step x asset arrays, variance reduction (antithetic variates,
   control variates using the closed-form Gibson-Schwartz forward price), and convergence
   diagnostics.

## Scope: what this deliberately will NOT attempt

Being explicit about this is itself part of the exercise — a real desk system differs from this
project in ways that are structural, not just a matter of more time:

- **No real-time/production market data infrastructure.** Forward curves, vol surfaces, and CDS
  curves will be either historical snapshots or parametrised by hand, not fed from a live pipeline.
- **No large-scale netting sets.** Real books can hold hundreds of heterogeneous trades with break
  clauses and bespoke CSA terms. This project uses a handful of trades to demonstrate the mechanics,
  not the scale.
- **No production-grade compute.** Real XVA engines run millions of paths on GPU/distributed
  clusters with incremental (trade-level delta) recomputation. This runs on a laptop at a much
  smaller path count, with variance reduction used to partially compensate.
- **No governance/control layer.** Model validation sign-off, daily PnL explain, regulatory capital
  (SA-CVA/BA-CVA), and front-office reconciliation are organisational processes outside the scope of
  a single-author technical project.
- **No genuine counterparty fundamentals for wrong-way risk.** Real WWR modelling draws on a
  counterparty's actual business exposure (e.g. an oil producer's revenue sensitivity to price).
  This project will only *illustrate* the mechanism via an assumed hazard-rate/price correlation,
  not model a real counterparty.

## Suggested architecture (mirrors the layering used in the PDE pricer)

```
MarketModel/        Gibson-Schwartz SDE, multi-commodity correlation, calibration stubs
PathSimulator/       Discretization schemes, streaming path generation (memory-aware)
Instrument/          Forward, European option, American/swing option payoffs
Portfolio/           Netting set aggregation
ExposureEngine/       EE / PFE / EPE / ENE, LSM regression for early exercise
Collateral/           Simplified CSA logic (threshold, MTA, MPOR)
XVACalculator/        CVA, FVA, (optional) simplified WWR overlay
VarianceReduction/    Antithetic variates, control variates
```

## Build roadmap

| Phase | Deliverable |
|---|---|
| 1 | Single-commodity Gibson-Schwartz simulator + exposure profile for a simple forward portfolio |
| 2 | Add option instruments with LSM for American/swing exposure |
| 3 | Multi-commodity correlation + netting across trade types |
| 4 | CVA/FVA calculation + simplified wrong-way risk |
| 5 | Performance: multithreading/vectorization, variance reduction, convergence diagnostics |

## How to talk about this in interviews

The honest framing is the strong framing: this project reproduces the *mathematics and structural
logic* of a commodities XVA desk (exposure definitions, netting, Monte Carlo necessity under high
dimensionality, LSM for early exercise, CVA/FVA integrals) at a scale and infrastructure level
appropriate for a single-author portfolio piece — not an attempt to replicate production
infrastructure, data feeds, or compute scale. Being able to state clearly what was simplified and why
demonstrates understanding of the real system, which is the point of the exercise.
