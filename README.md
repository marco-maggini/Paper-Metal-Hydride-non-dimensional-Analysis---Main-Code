MH_PCM.m is the core code. It comprises of three main sections:
a) data input, where geometry and thermo-physical properties of the system are defined;
b) absorption, where the hydrogen uptake process is simulated; and
c) desorption, where the hydrogen release is simulated.

inputData.m caters for the geometrical and physical properties of the system (i.e., hydride properties, diameter and length of the canister, PCM properties, reaction
  properties etc.).

absorption.m caters for the numerical simulation of the absorption process. It contains the set of 3 equations (1 ODE and 2 PDE's) which characterize the chemical reaction
  and the heat exchange within the canister. The set of algebraic equations is solved through ode15s.

desorptionWithStefanProblem.m caters for the numerical simulation of the desorption process. It contains the same 3 equations as of absorption, with slight differences.
It is solved through ode15s.

Secondary codes:
  getResults.m orders the data obtained from absorption.m and desorption.m into structs and substructs for practical reasons.

  absorptionCondition.m and desorptionStefanCondition.m contain the conditions for the 2 numerical simulations to stop (e.g., 99.5% of uptake complete).

  BaseDesign.mat contains a preset of 1-kWh canisters layouts with different L/D ratios.


FInally, "Sensitivity Analysis.zip" contains the .tikz of the sensitivity analysis.
