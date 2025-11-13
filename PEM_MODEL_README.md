# Advanced PEM Water Electrolyzer Model

A comprehensive computational model for Proton Exchange Membrane (PEM) water electrolyzers based on Franz et al. (2023) with full implementation of advanced features.

## Features

### ✓ Short-Term Enhancements
1. **Adaptive Mesh Refinement** - Higher resolution near catalyst layer interfaces
2. **Jacobian-Free Newton-Krylov Solver** - Advanced nonlinear solver for large systems
3. **Temperature-Dependent Water Content** - Springer correlation with T-dependent uptake

### ✓ Medium-Term Enhancements
4. **Two-Phase Flow in Porous Media** - Brooks-Corey model with Leverett J-function
5. **Dynamic Operation and Transients** - Time-dependent simulations
6. **Degradation Mechanisms** - Catalyst dissolution, membrane thinning, crossover increase

## Physics Included

- **Charge Transport**: Electronic and protonic conduction with Ohm's law
- **Mass Transport**: Multi-species diffusion (H2, O2, H2O) with concentration gradients
- **Electrochemical Kinetics**: Tafel kinetics for OER and HER
- **Water Management**: Water content evolution in ionomer with electro-osmotic drag
- **Two-Phase Flow**: Liquid water saturation with capillary effects
- **Degradation**: Time-dependent performance loss mechanisms

## Installation

### Requirements
```bash
pip install numpy scipy matplotlib
```

### Files
- `pem_electrolyzer_complete.py` - Main model implementation
- `test_pem_features.py` - Comprehensive test suite
- `CORRECTIONS_SUMMARY.md` - Detailed corrections documentation

## Quick Start

### Basic Usage

```python
from pem_electrolyzer_complete import Parameters, DynamicPEMElectrolyzer
import matplotlib.pyplot as plt

# 1. Create model with default parameters
params = Parameters()
model = DynamicPEMElectrolyzer(params)

# 2. Solve at single operating point
V_cell = 1.7  # Cell voltage [V]
result = model.solve_steady_state(V_cell, method='hybr')

if result['converged']:
    j = result['current_density']  # A/cm²
    print(f"Current density: {j:.4f} A/cm²")
```

## Model Parameters

### Operating Conditions
```python
params.T = 353.0           # Temperature [K] (80°C)
params.P_anode = 1e5       # Anode pressure [Pa] (1 bar)
params.P_cathode = 1e5     # Cathode pressure [Pa]
```

### Membrane Properties
```python
params.membrane_type = 'N117'
params.EW = 1100.0         # Equivalent weight [g/mol]
params.λ = 14.0            # Water content [-]
params.L_mem_dry = 178e-6  # Dry thickness [m]
```

### Catalyst Layers
```python
# Anode (OER)
params.L_aCL = 6e-6        # Thickness [m]
params.loading_aCL = 1.4   # Catalyst loading [mg/cm²]
params.SA_IrOx = 100.0     # Specific area [m²/g]
params.j0_aCL_ref = 1e-7   # Exchange current density [A/m²]
params.b_aCL_ref = 0.042   # Tafel slope [V]

# Cathode (HER)
params.L_cCL = 12e-6
params.loading_cCL = 0.8
params.SA_PtC = 400.0
params.j0_cCL_ref = 1.0
params.b_cCL_ref = 0.030
```

## Testing

Run comprehensive test suite:
```bash
python test_pem_features.py
```

## References

1. **Franz, A., et al. (2023)**. "Continuum modeling of PEM water electrolysis: A comprehensive review and comparison of approaches." *Journal of Power Sources*, 559, 232582.

2. **Springer, T.E., et al. (1991)**. "Polymer electrolyte fuel cell model." *Journal of the Electrochemical Society*, 138(8), 2334-2342.

For detailed corrections and implementation notes, see `CORRECTIONS_SUMMARY.md`.
