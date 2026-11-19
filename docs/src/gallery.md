# Gallery: resilience across ecorecipes

ResilienceCore computes metrics on operators. The modelling packages produce those operators
from their own models, and each one shows how in its own documentation. This page links to them,
grouped by the kind of object the metrics are computed on. The table on the
[home page](index.md) says which metrics apply to which object.

## Simulation and SciML models

[DynamicalResilience.jl](https://ecorecipes.github.io/DynamicalResilience.jl/dev/) brings
ResilienceCore to SciML problems:
- steady states, R₀, Floquet multipliers and continuation to a fold or Hopf point;
- simulation-based recovery, basin stability and early-warning indicators, on SciML's integrator
  interface.

These packages use it:

| Package | Tutorial | What it shows |
|---|---|---|
| Chytridiomycosis.jl | [Thresholds and resilience](https://ecorecipes.github.io/Chytridiomycosis.jl/dev/tutorials/07_resilience/) | resilience falls to zero at the Hopf threshold, and a load-tipping fold |
| PartiallySpecifiedModels.jl | [Resilience of fitted models](https://ecorecipes.github.io/PartiallySpecifiedModels.jl/dev/tutorials/42_resilience/) | uncertainty in resilience and the Hopf point from bootstrap fits of a mechanistic and a spline model |
| PhysiologicallyBasedDemographicModels.jl | [Temperature, annual cycles and recovery](https://ecorecipes.github.io/PhysiologicallyBasedDemographicModels.jl/dev/tutorials/65_resilience/) | resilience against temperature, annual maps and recovery under real weather |
| ArkModelZoo.jl | [Resilience of agent-based models](https://ecorecipes.github.io/ArkModelZoo.jl/dev/tutorials/resilience/resilience/) | Daisyworld's basin and early warnings before collapse; Wolf–Sheep recovery and its noise floor |
| EcologicalCellularAutomata.jl | [Resilience of a viral front near its threshold](https://ecorecipes.github.io/EcologicalCellularAutomata.jl/dev/tutorials/34_viral_resilience/) | return time and establishment near a critical point, and warnings along a ramp |
| IndividualBasedPopulationDynamics.jl | [Resilience of an individual-based population](https://ecorecipes.github.io/IndividualBasedPopulationDynamics.jl/dev/tutorials/04_resilience/) | return of a size distribution after a cull, against its projection model and the demographic-noise floor |

## Equilibria of community models

| Package | Tutorial | What it shows |
|---|---|---|
| CommunityMatrices.jl | [Plant communities](https://ecorecipes.github.io/CommunityMatrices.jl/dev/tutorials/01_plant_communities/) | structural stability and sign patterns of empirical interaction matrices against shuffled nulls |
| CommunityMatrices.jl | [Gut microbiota](https://ecorecipes.github.io/CommunityMatrices.jl/dev/tutorials/02_gut_microbiota/) | a fitted Lotka–Volterra model: assembly, invasion, an antibiotic pulse and S-map Jacobians |
| MatrixCommunityModels.jl | [Beetles and riparian guilds](https://ecorecipes.github.io/MatrixCommunityModels.jl/dev/tutorials/26_resilience/) | stage-structured communities, periodic environments, and two ways of pressing mortality |
| DominguezGarcia2019.jl | [Tutorials](https://ecorecipes.github.io/DominguezGarcia2019.jl/dev/) | a reproduction of Domínguez-García, Dakos & Kéfi (2019): 27 metrics on simulated and observed food webs, and which of them measure the same thing |

## Structured populations

- **MatrixProjectionModels.jl and IntegralProjectionModels.jl:** `linearise` comes from
  extensions in ResilienceCore. It gives a `ProjectionMatrix`, or a `Generator` with mesh
  weights.
- **FiniteStatePopulationDynamics.jl:** generator problems give a `Generator`, the Jacobian at an
  affine equilibrium, or a `PeriodicOperator`
  ([section](https://ecorecipes.github.io/FiniteStatePopulationDynamics.jl/dev/#Resilience-operators)).
- **ContinuousStatePopulationDynamics.jl:** the same, with mesh weights
  ([section](https://ecorecipes.github.io/ContinuousStatePopulationDynamics.jl/dev/#Resilience-operators)).

## Markov chains and sign patterns

| Package | Tutorial | What it shows |
|---|---|---|
| FiniteKernels.jl | [Markov chains and their resilience](https://ecorecipes.github.io/FiniteKernels.jl/dev/tutorials/02_markov_chains/) | `MarkovKernel` metrics on Horn's forest succession |
| CausalLoopDiagrams.jl | [Qualitative stability](https://ecorecipes.github.io/CausalLoopDiagrams.jl/dev/tutorials/02_qualitative_stability/) | what the signs of a diagram alone say about stability and press responses |
