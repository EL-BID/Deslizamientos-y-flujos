# Input Parameters Reference

<img width="1755" height="1007" alt="image" src="https://github.com/user-attachments/assets/9cdc9046-57fe-48ce-81c8-292e16388103" />

## Section A: Model Identification

### Name of Data File
Base filename for all output files.

- The plugin automatically appends the `.DAT` extension if not provided.
- The filename should follow standard naming conventions and avoid special characters.

### Description
Case study identifier that appears in the output file header.

**Recommended format:**

```text
Location-Type-Year
```

**Example:**

```text
Vajont-Landslide-1963
```

### Modeling Types

#### Deterministic Propagation Model
- Uses fixed parameters.
- Performs a single simulation run.
- Suitable for well-characterized scenarios.

#### Probabilistic Propagation Model
- Incorporates parameter uncertainty.
- Restricted in the public version.

#### Probabilistic Initiation Model
- Focuses on failure probability.
- Restricted in the public version.

### Landslide Types

#### Mudflow (One-phase)
- Single-phase solid-dominated flow.
- Suitable for debris flows with high sediment concentration.

#### Flood (One-phase)
- Single-phase fluid flow.
- Suitable for clear-water flooding scenarios.
- Restricted in the public version.

#### Flowslide (One-phase + Effective Density)
- One-phase model with effective density accounting for buoyancy.
- Restricted in the public version.

#### Debris Flow (Two-phase)
- Simulates solid-fluid interaction.
- Suitable for sediment-laden flows.

#### Landslide in Reservoir (Two-phase)
- Two-phase submerged flow model.
- Suitable for submarine landslides.
- Restricted in the public version.

---

## Section B: Source and Boundary Conditions

### Random Factor
Stochastic factor used for probabilistic simulations.

- Range: `0 – 1`
- Default: `0` (deterministic simulation)

### Injection Types

#### No Injection
- Standard initial condition without mass addition.

#### Normal Injection
- Mass addition at a specified rate.
- Restricted.

#### Riemann Injection
- Characteristic-based inflow condition.
- Restricted.

### Phases

#### Solid
- Granular material phase.
- Enabled depending on the selected landslide type.

#### Fluid
- Water or fluid phase.
- Enabled depending on the selected landslide type.

### Inserting Walls
Boundary condition used to define domain walls.

> Partially functional.

### Source Points

#### Source Points (3D)
- Reads particle positions from an external `.PTS` file.

#### Dam-break (2D)
- Rectangular reservoir configuration.
- Restricted.

#### Cloud (3D)
- Random particle distribution.
- Restricted.

#### Cylinder (3D)
- Cylindrical particle arrangement.
- Functional.

#### Injection (3D)
- Continuous particle injection.
- Restricted.

#### Ellipsoid (3D)
- Ellipsoidal particle cloud.
- Restricted.

#### Restart (3D)
- Continues a simulation from a previous run.
- Restricted.

#### 2D to 3D
- Extrudes a 2D configuration into 3D.
- Restricted.

---

## Section C: Numerical Methods

### Particle Approximation
Method used to calculate gradients in the SPH equations.

### NNPS
Nearest Neighbor Particle Search algorithm.

### Smoothing Length
Method used to update the kernel support radius.

### SKF
Kernel function used for interpolation weighting.

### Density Summation
Method used to calculate the density field.

### Average Velocity
Option for velocity field smoothing.

### Normalized Density
Density normalization approach used to improve numerical stability.

---

## Section D: Physical Parameters

### Gravitational Acceleration
Standard gravity constant.

- Units: `m/s²`
- Default: `9.81`
- Typical range: `9.76 – 9.83`

### Mixture Density
Bulk density of the flow mixture.

- Units: `kg/m³`

### Solid Density
Density of solid particles.

- Units: `kg/m³`
- Typical range: `2400 – 2700`

### Fluid Density
Density of the fluid phase.

- Units: `kg/m³`
- Typical value for water: `1000`

### Manning's Roughness Coefficient
Bottom friction parameter.

- Units: `s/m¹/³`
- Range: `0.01 – 0.15`
  - `0.01` = very smooth
  - `0.15` = very rough

### Voellmy Turbulent Coefficient
Turbulent friction parameter.

- Units: `m/s²`
- Range: `500 – 5000`

### Cohesion
Material cohesion for the Bingham rheological model.

- Units: `Pa`
- Range: `0 – 50,000`

### Shear Strength
Yield stress for the Bingham model.

- Units: `Pa`
- Range: `0 – 100,000`

### Viscosity
Dynamic viscosity for Newtonian or Bingham models.

- Units: `Pa·s`

Typical values:

| Material | Viscosity (Pa·s) |
|-----------|----------------|
| Water | 0.001 |
| Debris flow | 0.1 - 10 |

### Friction Angle
Basal friction angle used in frictional models.

- Entered as the tangent value.
- Range:
  - `tan(0°) = 0`
  - `tan(45°) = 1`

### Minimum Flow Thickness
Lower limit for flow depth to prevent numerical instability.

- Units: `m`
- Default: `0.01`

### Initial Relative Height
Initial flow height in the source region.

- Units: `m`

### Volumetric Stiffness
Bulk modulus for two-phase flows.

- Units: `N/m²`
- Range: `10⁶ – 10⁹`

### Porosity Lower Limit
Minimum porosity allowed to avoid unphysical values.

- Default: `0.01`

---

## Section E: Advanced Parameters

### Erosion Law
Selection of the erosion and entrainment model.

### Drag Law
Fluid-solid interaction model for two-phase simulations.

### Coordinate System
Defines whether calculations use global or local coordinates.

### Coarse Mesh
Multi-resolution modelling capability.

> Restricted in the public version.

### Check Points
Monitoring point definition for simulation tracking.

> Restricted in the public version.

### Output Variables
GID mask selection used to define result file outputs.
