# Metric catalogue

Every resilience metric in ResilienceCore and its sister packages, by name. The table on the
[home page](index.md) says which operators each metric applies to.

## The 27 metrics of Domínguez-García, Dakos & Kéfi (2019)

[`stability_profile`](@ref) computes all 27 for any community that implements the community
interface. [`jacobian_profile`](@ref) gives the 11 that need only a Jacobian. Values keep their
natural orientation; the paper flips signs so that larger always means more stable. The last
column gives the paper's value in terms of ours (N is the number of species present).

| Paper | Field | Metric | ResilienceCore | Paper value |
|---|---|---|---|---|
| R0 | `R0` | reactivity | [`reactivity`](@ref) | −R0 |
| Is | `Is` | intrinsic stochastic invariability | [`stochastic_invariability`](@ref) | Is |
| Amax | `Amax` | maximal amplification | [`max_amplification`](@ref) | −Amax |
| tmax | `tmax` | time of maximal amplification | [`max_amplification`](@ref) | −tmax |
| R∞ | `Rinf` | asymptotic resilience | [`resilience`](@ref) | R∞ |
| MR0 | `MR0` | median reactivity | [`median_reactivity`](@ref) | −MR0 |
| MAmax | `MAmax` | median maximal amplification | [`median_max_amplification`](@ref) | −MAmax |
| Mtmax | `Mtmax` | time of median maximal amplification | [`median_max_amplification`](@ref) | −Mtmax |
| MR∞ | `MRinf` | median resilience | [`median_resilience`](@ref) | MR∞ |
| S | `S` | community sensitivity `Σᵢ\|Σⱼ sᵢⱼ\|` | [`press_sensitivity`](@ref) | −S |
| ⟨sᵢⱼ⟩ | `mean_sij` | mean absolute sensitivity | [`total_sensitivity`](@ref) / N² | −⟨sᵢⱼ⟩ |
| TM^G | `TM` | global mortality tolerance | [`global_mortality_tolerance`](@ref) | TM |
| ⟨TM^L⟩ | `mean_TM` | mean local mortality tolerance | [`local_mortality_tolerances`](@ref) | ⟨TM⟩ |
| TM^L_min | `min_TM` | least local mortality tolerance | [`local_mortality_tolerances`](@ref) | min TM |
| ⟨CE⟩ | `mean_CE` | mean secondary extinctions per removal | [`removal_responses`](@ref) | −N/(N−1) ⟨CE⟩ |
| ⟨TE⟩ | `TE` | fraction removed before half are lost | [`extinction_tolerance`](@ref) / N | TE |
| CE_max | `max_CE` | most secondary extinctions per removal | [`removal_responses`](@ref) | −N/(N−1) max CE |
| RM^G | `RM` | resistance to a global press | [`global_press_response`](@ref) | −RM |
| ⟨RM^L⟩ | `mean_RM` | resistance to local presses, mean | [`local_press_responses`](@ref) | −⟨RM⟩ |
| RM^L_max | `max_RM` | resistance to local presses, worst | [`local_press_responses`](@ref) | −max RM |
| SM^G | `SM` | L1 sensitivity to a global press | [`global_press_response`](@ref) | −SM |
| ⟨SM^L⟩ | `mean_SM` | L1 sensitivity to local presses, mean | [`local_press_responses`](@ref) | −⟨SM⟩ |
| SM^L_max | `max_SM` | L1 sensitivity to local presses, worst | [`local_press_responses`](@ref) | −max SM |
| ⟨RE⟩ | `mean_RE` | resistance to one removal, mean | [`removal_responses`](@ref) | −⟨RE⟩ |
| RE_max | `max_RE` | resistance to one removal, worst | [`removal_responses`](@ref) | −max RE |
| ⟨SE⟩ | `mean_SE` | L1 sensitivity to one removal, mean | [`removal_responses`](@ref) | −⟨SE⟩ |
| SE_max | `max_SE` | L1 sensitivity to one removal, worst | [`removal_responses`](@ref) | −max SE |

The published Dataset S1 differs from the SI's definitions in three places: the norm of
`Amax`/`tmax`, the pulse weights of the median metrics, and the size of the global press.
DominguezGarcia2019.jl reproduces those conventions.

## Other metrics in ResilienceCore

| Metric | Function | Objects |
|---|---|---|
| asymptotic growth rate λ | [`growth_rate`](@ref) | projections, periodic projections |
| damping ratio | [`resilience`](@ref) | projections, generators, Markov kernels |
| stability | [`isstable`](@ref) | Jacobians |
| transient envelope | [`amplification_envelope`](@ref) | Jacobians, projections |
| sensitivity and elasticity of λ | [`growth_sensitivity`](@ref), [`elasticity`](@ref) | projections, generators |
| stationary distribution, mixing time | [`stationary_distribution`](@ref), [`mixing_time`](@ref) | Markov kernels |
| qualitative (sign) stability | [`quasi_sign_stability`](@ref) | sign patterns |
| sign determinacy of press responses | [`press_sign_determinacy`](@ref) | sign patterns |

## Metrics in other packages

| Metric | Function | Package |
|---|---|---|
| basic reproduction number R₀ | `reproduction_number` | DynamicalResilience.jl |
| Floquet multipliers of a cycle | `period_map_jacobian`, `periodic_orbit` | DynamicalResilience.jl |
| resilience along a parameter, to a fold or Hopf point | `resilience_along` | DynamicalResilience.jl |
| simulated recovery time | `recovery` | DynamicalResilience.jl |
| basin stability | `basin_stability` | DynamicalResilience.jl |
| early-warning indicators (variance, lag-1 autocorrelation, Kendall τ) | `early_warning` | DynamicalResilience.jl |
| structural stability (feasibility domain and angle) | `feasibility_domain`, `feasibility_angle` | CommunityMatrices.jl |
| empirical invariability, S-map Jacobians | `invariability`, `smap_jacobians` | CommunityMatrices.jl |
