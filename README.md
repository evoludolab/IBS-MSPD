# IBS-MSPD

Mutualistic spatial prisoner's dilemma games

Individual-based simulations for modelling the rich and complex dynamics of cooperation between two species engaging in the spatial donation game. Each species is arranged on its own lattice layer. Interactions occur between species and hence across layers, while competition arises within species and layers.

This source code reproduces the simulation results published in Hauert, C & Szabó, G. (2024) *Spontaneous symmetry breaking of cooperation between species* PNAS Nexus, 3 (9) pgae326, https://doi.org/10.1093/pnasnexus/pgae326

## Source code

#### ch2layer.c
Donation game on a two-layer square lattice. Interactions occur with neighbours on the *other* layer. Stochastic imitation of a random neighbor happens within the *same* layer. The Monte-Carlo (MC) simulations investigate a sequence of critical phase transitions for decreasing cost-to-benefit ratios *r*:
1. for large *r* cooperation is too costly and defection dominates; 
2. when lowering *r* cooperators survive at equal frequencies in both species. This transition is accompanied by diverging correlations and fluctuations;
3. lowering *r* further results in intriguing, spontaneous symmetry breaking of cooperation between species. This transition is reminiscent of spontaneous magnetization in ferromagnets; 
4. finally, for small *r*, bursts of defection cause diverging fluctuations, which drive the populations into absorbing states with perfect asymmetry with one species cooperating and the other defecting. 
Returns the frequency of cooperators in each layer and their fluctuations together with the average total population payoff as a function of the cost-to-benefit ratio, *r*.

#### ch2layrp.c
Same as `ch2layer.c` but for starting from a prepared initial state.

#### ch2lcov.c</dt>
Same as `ch2layer.c` but additionally returns the average payoff for each population as well as the covariance of cooperation between the two layers.

#### chss.c</dt>
Same as `ch2layer.c` but for studying the pattern evolution in a small box, starting from a random initial state. Snapshots are saves as an eps file.

#### opfluct.c</dt>
Same as `ch2layer.c` but additionally returns the order parameter and its fluctuations as well as the covariance of cooperation between the two layers.


## Simulation data

#### ch2layer.dat
Baseline simulation data from `ch2layer.c` on `1000x1000` lattice.

#### ch2layer.dxt
Aggregated data from several runs of `ch2layer.c` with different parameters (longer relaxation and sampling times as well as lattice sizes up to `1800x1800`).

#### ch2layerp.dxt
Aggregated data from several runs of `ch2layerp.c` with different parameters (longer relaxation and sampling times as well as lattice sizes up to `1000x1000`).

#### ch2lcov.dat
Baseline simulation data from `ch2lcov.c` on `300x300` lattice.

#### ch2lcov.dxt
Aggregated data from several runs of `ch2lcov.c` with different parameters (longer relaxation and sampling times as well as lattice sizes up to `800x800`).

#### chss022.eps
Sample output from `chss.c` for a cost/benefit ratio of `r=0.0022`.

#### opfluct.dat</dt>
Baseline simulation data from `opfluct.c` on `800x800` lattice.

#### opfluct.dxt</dt>
Aggregated data from several runs of `opfluct.c` with different parameters (longer relaxation and sampling times as well as lattice sizes up to `1800x1800`).


## EvoLudo project

For interactive online tutorials and visualizations of the different dynamical regimes check out https://www.evoludo.org. 
