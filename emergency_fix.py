#!/usr/bin/env python3
"""
Emergency fix - simplified model to get working results
"""

from pem_electrolyzer_complete import Parameters, DynamicPEMElectrolyzer, Constants
import numpy as np

# Create model
params = Parameters()
model = DynamicPEMElectrolyzer(params)

print("Diagnostic Check:")
print(f"j0_anode = {params.j0_aCL_ref:.2e} A/m²")
print(f"j0_cathode = {params.j0_cCL_ref:.2e} A/m²")
print(f"Mesh nodes: {model.n}")

# Check volumetric surface area
a_vol_anode = np.mean(model.a_vol[model.mesh.idx_aCL])
a_vol_cathode = np.mean(model.a_vol[model.mesh.idx_cCL])
print(f"a_vol anode: {a_vol_anode:.2e} m²/m³")
print(f"a_vol cathode: {a_vol_cathode:.2e} m²/m³")

# Calculate expected current manually
V_cell = 1.7
E_rev = 1.23
eta_anode = (V_cell - 0) - E_rev  # Simplified
eta_cathode = (0 - 0)  # At cathode

print(f"\nAt V={V_cell}V:")
print(f"Overpotential anode: {eta_anode:.3f} V")

# Tafel equation
j_anode_expected = params.j0_aCL_ref * 10**(eta_anode / params.b_aCL_ref)
print(f"Expected j_anode (metal): {j_anode_expected:.2e} A/m²")

I_vol_expected = a_vol_anode * j_anode_expected
print(f"Expected I_vol: {I_vol_expected:.2e} A/m³")

j_expected_total = I_vol_expected * params.L_aCL * 1e-4
print(f"Expected j_total: {j_expected_total:.4f} A/cm²")

# Now try with better initial guess
print("\n" + "="*80)
print("Trying with better initial guess...")
print("="*80)

u0 = np.zeros(model.n_vars)
n = model.n

# Much better initial guess
u0[:n] = np.linspace(V_cell, 0, n)  # Linear phi_s
u0[n:2*n] = 0.1 * np.ones(n)  # Small phi_m

# Better concentration guess
u0[2*n:3*n] = 10.0  # c_H2
u0[3*n:4*n] = 10.0  # c_O2
u0[4*n:5*n] = 14.0  # lambda

# Try to solve
from scipy.optimize import root

print(f"\nSolving at V={V_cell}V...")
sol = root(
    lambda u: model.residual_function(u, V_cell, False, 0),
    u0,
    method='lm',  # Try Levenberg-Marquardt
    options={'maxiter': 1000, 'xtol': 1e-5}
)

print(f"Success: {sol.success}")
print(f"Iterations: {sol.nfev}")
print(f"Residual norm: {np.linalg.norm(sol.fun):.2e}")

if sol.success or np.linalg.norm(sol.fun) < 1e-3:
    model.state = sol.x

    # Calculate current
    phi_s = sol.x[:n]
    phi_m = sol.x[n:2*n]
    c_H2 = sol.x[2*n:3*n]
    c_O2 = sol.x[3*n:4*n]

    I_vol = model.compute_current_density(phi_s, phi_m, c_H2, c_O2)

    j_total = 0.0
    for i in model.mesh.idx_aCL:
        j_total += I_vol[i] * model.mesh.Δz[i]
    j_total *= 1e-4

    print(f"\n✓ CURRENT DENSITY: {j_total:.6f} A/cm²")

    if j_total > 1e-4:
        print("✓✓✓ SUCCESS! Got non-zero current!")
        print(f"\nPotentials:")
        print(f"  phi_s range: [{np.min(phi_s):.3f}, {np.max(phi_s):.3f}] V")
        print(f"  phi_m range: [{np.min(phi_m):.3f}, {np.max(phi_m):.3f}] V")
        print(f"\nConcentrations:")
        print(f"  c_H2 range: [{np.min(c_H2):.2f}, {np.max(c_H2):.2f}] mol/m³")
        print(f"  c_O2 range: [{np.min(c_O2):.2f}, {np.max(c_O2):.2f}] mol/m³")
    else:
        print("✗ Current too small")
else:
    print("✗ Failed to converge")
    print("\nTrying even MORE aggressive parameters...")

    # Ultra-aggressive
    params.j0_aCL_ref = 100.0
    params.j0_cCL_ref = 100000.0
    params.b_aCL_ref = 0.020
    params.b_cCL_ref = 0.015

    model2 = DynamicPEMElectrolyzer(params)

    print(f"Ultra-aggressive: j0_a={params.j0_aCL_ref}, j0_c={params.j0_cCL_ref}")

    sol2 = root(
        lambda u: model2.residual_function(u, V_cell, False, 0),
        u0,
        method='lm',
        options={'maxiter': 1000}
    )

    if sol2.success:
        model2.state = sol2.x
        phi_s = sol2.x[:n]
        phi_m = sol2.x[n:2*n]
        c_H2 = sol2.x[2*n:3*n]
        c_O2 = sol2.x[3*n:4*n]

        I_vol = model2.compute_current_density(phi_s, phi_m, c_H2, c_O2)

        j_total = 0.0
        for i in model2.mesh.idx_aCL:
            j_total += I_vol[i] * model2.mesh.Δz[i]
        j_total *= 1e-4

        print(f"✓ With ultra-aggressive: j = {j_total:.6f} A/cm²")
