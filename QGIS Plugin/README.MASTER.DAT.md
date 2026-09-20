# Simulation Parameters Reference (MASTER.DAT)

<img width="1757" height="1008" alt="image" src="https://github.com/user-attachments/assets/dd2f19e2-94a6-4fe7-8737-803b3b8a858a" />

## Section A: Simulation Control

### Name of Data File
Base filename matching the corresponding `.DAT` file.

### Description
Simulation identifier included in the master file header.

### Modeling Types
Same modeling options available in **Tab 1**, used for multi-run coordination and simulation setup.

### Problem Types

#### Shallow-Water Equation
- Uses a depth-averaged flow formulation.
- Suitable for large-scale flow simulations where vertical velocity variations are negligible.

#### Fractional-Step Method
- Employs an operator-splitting approach.
- Solves complex governing equations in separate computational steps.

### Algorithm Types

#### 4th Order Runge-Kutta
- High-order explicit time integration method.
- Provides high accuracy but requires restrictive time steps for stability.

#### Semi-Implicit
- Combines explicit and implicit integration schemes.
- Suitable for stiff problems and allows larger time steps.

#### Slow Flow
- Optimized for low-velocity flow conditions.

#### Two-Step Taylor-Galerkin
- Second-order accurate numerical scheme.
- Well suited for wave propagation and transient flow problems.

#### Random Walk
- Stochastic particle-tracking approach.
- Used for uncertainty and dispersion analyses.

#### Monte Carlo Method
- Statistical sampling method used for uncertainty quantification.
- Requires multiple simulation realizations.

#### First-Order Second-Moment Method
- Analytical uncertainty propagation technique.
- Estimates output variability using statistical moments.

---

## Section B: Time Parameters

### Time-Step
Initial simulation time-step size.

- Units: `s`
- Should satisfy the CFL stability condition.

### Total Analysis Time
Total duration of the simulation.

- Units: `s`

### Maximum Number of Time-Steps
Maximum number of iterations allowed.

- Prevents simulations from running indefinitely.

### Print at Every
Frequency of console output.

- Units: simulation steps

### Save at Every
Frequency of result file generation.

- Units: simulation steps

### Plot at Every
Frequency of visualization updates.

- Units: simulation steps

### Time-Stepping Technique

#### Constant
- Uses a fixed time-step throughout the simulation.

#### Adaptive
- Automatically adjusts the time-step according to stability and accuracy criteria.

---

## Section C: Probabilistic Parameters

### Uncertain Variables

Parameters that can be treated as random variables during uncertainty analyses.

#### Mixture Density
- Represents uncertainty in flow density.

#### Cohesion
- Represents variability in material cohesion.

#### Voellmy Turbulent Coefficient
- Represents uncertainty in turbulent friction resistance.

#### Viscosity
- Represents uncertainty in flow rheology and resistance.

#### Basal Friction Angle
- Represents uncertainty in slope friction and flow mobility.

### Mean
Expected value assigned to each uncertain variable.

### Standard Deviation
Measure of variability around the mean value.

### Probability Distribution

#### Triangular Distribution
- Suitable for bounded parameters with a known minimum, maximum, and most likely value.

#### Normal (Gaussian) Distribution
- Suitable for symmetric uncertainty around the mean.

#### Log-Normal Distribution
- Suitable for positively skewed parameters that cannot take negative values.

### Random Number
Seed value for the random number generator.

- Ensures reproducibility of probabilistic simulations.

---

## Section D: PGA and Relative Height (Conditional Options)

### Activate PGA
Enables probabilistic analysis based on ground-motion intensity.

> Restricted in the public version.

### Activate Relative Height
Enables simulations that account for water-level variability.

> Restricted in the public version.

### Return Period Values
Statistical recurrence intervals used in probabilistic analyses.

> Restricted in the public version.

### Number of Years
Time horizon used for probability calculations.

> Restricted in the public version.
