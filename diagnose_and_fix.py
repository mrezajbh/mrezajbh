#!/usr/bin/env python3
"""
Diagnostic and aggressive parameter tuning for PEM electrolyzer
"""

from pem_electrolyzer_complete import Parameters, DynamicPEMElectrolyzer
import numpy as np
import matplotlib.pyplot as plt

def check_model_properties(model):
    """Diagnose model properties"""
    print("\n" + "=" * 80)
    print("MODEL DIAGNOSTICS")
    print("=" * 80)

    # Check volumetric surface area
    if len(model.mesh.idx_aCL) > 0:
        a_vol_anode = np.mean(model.a_vol[model.mesh.idx_aCL])
        print(f"Anode volumetric area: {a_vol_anode:.2e} m²/m³")
        if a_vol_anode < 1e6:
            print("  ⚠ WARNING: Very low! Should be > 1e6 m²/m³")
    else:
        print("  ✗ No anode CL nodes!")

    if len(model.mesh.idx_cCL) > 0:
        a_vol_cathode = np.mean(model.a_vol[model.mesh.idx_cCL])
        print(f"Cathode volumetric area: {a_vol_cathode:.2e} m²/m³")
        if a_vol_cathode < 1e6:
            print("  ⚠ WARNING: Very low! Should be > 1e6 m²/m³")
    else:
        print("  ✗ No cathode CL nodes!")

    # Check conductivities
    if len(model.mesh.idx_mem) > 0:
        sigma_m = np.mean(model.σ_m[model.mesh.idx_mem])
        print(f"Membrane conductivity: {sigma_m:.2f} S/m")
        if sigma_m < 1.0:
            print("  ⚠ WARNING: Low conductivity!")
    else:
        print("  ✗ No membrane nodes!")

    # Check layer thicknesses
    print(f"\nLayer thicknesses:")
    print(f"  PTL: {model.params.L_PTL*1e6:.1f} μm ({len(model.mesh.idx_PTL)} nodes)")
    print(f"  Anode CL: {model.params.L_aCL*1e6:.2f} μm ({len(model.mesh.idx_aCL)} nodes)")
    print(f"  Membrane: {model.params.L_mem_wet*1e6:.1f} μm ({len(model.mesh.idx_mem)} nodes)")
    print(f"  Cathode CL: {model.params.L_cCL*1e6:.2f} μm ({len(model.mesh.idx_cCL)} nodes)")
    print(f"  GDL: {model.params.L_GDL*1e6:.1f} μm ({len(model.mesh.idx_GDL)} nodes)")

def get_very_aggressive_parameters():
    """Return very aggressive parameters for testing"""
    params = Parameters()

    # Operating conditions
    params.T = 353.0  # 80°C
    params.P_anode = 1e5
    params.P_cathode = 1e5

    # VERY AGGRESSIVE: Much higher exchange currents
    params.j0_aCL_ref = 0.1  # 0.1 A/m² (was 1e-7)
    params.j0_cCL_ref = 10000.0  # 10,000 A/m² (was 1.0)

    # Lower Tafel slopes = faster kinetics
    params.b_aCL_ref = 0.025  # V
    params.b_cCL_ref = 0.020  # V

    # Much higher catalyst loadings
    params.loading_aCL = 4.0  # mg/cm²
    params.loading_cCL = 2.0  # mg/cm²

    # Higher specific areas
    params.SA_IrOx = 200.0  # m²/g
    params.SA_PtC = 800.0   # m²/g

    # Higher electronic conductivities
    params.σ_e_aCL = 3000.0  # S/m
    params.σ_e_cCL = 600.0   # S/m

    return params

def simple_test():
    """Simple test with very aggressive parameters"""
    print("=" * 80)
    print("AGGRESSIVE PARAMETER TEST")
    print("=" * 80)

    params = get_very_aggressive_parameters()

    print(f"\nVery Aggressive Parameters:")
    print(f"  j0_anode:    {params.j0_aCL_ref:.2e} A/m²")
    print(f"  j0_cathode:  {params.j0_cCL_ref:.2e} A/m²")
    print(f"  b_anode:     {params.b_aCL_ref:.3f} V")
    print(f"  b_cathode:   {params.b_cCL_ref:.3f} V")
    print(f"  Loading (a): {params.loading_aCL:.1f} mg/cm²")
    print(f"  Loading (c): {params.loading_cCL:.1f} mg/cm²")
    print(f"  SA IrOx:     {params.SA_IrOx:.0f} m²/g")
    print(f"  SA PtC:      {params.SA_PtC:.0f} m²/g")

    # Create model
    model = DynamicPEMElectrolyzer(params)

    # Run diagnostics
    check_model_properties(model)

    # Try solving
    print("\n" + "=" * 80)
    print("ATTEMPTING TO SOLVE")
    print("=" * 80)

    results = []
    for V in [1.5, 1.6, 1.7, 1.8, 1.9, 2.0]:
        print(f"\nTrying V = {V:.2f} V...")
        result = model.solve_steady_state(V, method='hybr')

        if result['converged']:
            j = result['current_density']
            print(f"  ✓ Converged: j = {j:.6f} A/cm²")

            if j > 1e-6:  # Non-zero current
                results.append((V, j))
                print(f"    → SUCCESS! Non-zero current achieved!")
            else:
                print(f"    ⚠ Zero current despite convergence")
        else:
            print(f"  ✗ Failed to converge")

    # Report results
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)

    if len(results) > 0:
        print(f"\n✓ SUCCESS! Got {len(results)} converged points with non-zero current:")
        for V, j in results:
            print(f"  V = {V:.2f} V → j = {j:.6f} A/cm²")

        # Plot
        V_vals = [r[0] for r in results]
        j_vals = [r[1] for r in results]

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(j_vals, V_vals, 'o-', linewidth=2, markersize=10)
        ax.set_xlabel('Current Density [A/cm²]', fontsize=12)
        ax.set_ylabel('Cell Voltage [V]', fontsize=12)
        ax.set_title('PEM Electrolyzer - Aggressive Parameters', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)
        plt.tight_layout()
        plt.savefig('aggressive_polarization.png', dpi=150)
        print(f"\n✓ Plot saved to: aggressive_polarization.png")

        # Calculate efficiency
        if max(j_vals) >= 1.0:
            V_at_1A = np.interp(1.0, j_vals, V_vals)
            eff = 1.23 / V_at_1A * 100
            print(f"\n✓ At 1 A/cm²:")
            print(f"  Voltage: {V_at_1A:.3f} V")
            print(f"  Efficiency: {eff:.1f}%")

    else:
        print("\n✗ STILL FAILING")
        print("\nPossible issues:")
        print("1. Exchange current densities still too low")
        print("2. Volumetric surface area calculation issue")
        print("3. Initial guess problem")
        print("4. Numerical precision issues")
        print("\nRecommendation: Increase j0_aCL_ref to 1.0 or even 10.0")

    return model, results

if __name__ == "__main__":
    model, results = simple_test()

    # Additional suggestions
    if len(results) == 0:
        print("\n" + "=" * 80)
        print("NEXT STEPS TO TRY")
        print("=" * 80)
        print("\n1. EXTREME PARAMETERS:")
        print("   params.j0_aCL_ref = 10.0  # Very high!")
        print("   params.j0_cCL_ref = 100000.0")
        print("   params.b_aCL_ref = 0.015")
        print("   params.b_cCL_ref = 0.015")

        print("\n2. SIMPLIFY THE MODEL:")
        print("   - Reduce mesh refinement")
        print("   - Use simpler initial guess")
        print("   - Try 'lm' method instead of 'hybr'")

        print("\n3. CHECK CALCULATION:")
        print("   # Manually calculate expected current")
        print("   η = 0.3  # 300 mV overpotential")
        print("   j_expected = j0 * 10^(η/b)")
        print("   print(f'Expected j: {j_expected:.2e} A/m²')")

        print("\n4. ALTERNATIVE APPROACH:")
        print("   Consider using literature parameters from")
        print("   specific experimental papers on PEM electrolyzers")
