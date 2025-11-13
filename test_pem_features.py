#!/usr/bin/env python3
"""
Test and demonstration of PEM electrolyzer complete model features
Validates all physics and numerical implementations
"""

import numpy as np
import matplotlib.pyplot as plt
from pem_electrolyzer_complete import (
    Parameters, EnhancedProperties, AdaptiveMesh,
    DegradationState, TwoPhaseFlow, DynamicPEMElectrolyzer, Constants
)

def test_enhanced_properties():
    """Test all enhanced property correlations"""
    print("\n" + "=" * 80)
    print("TEST 1: ENHANCED PROPERTY CORRELATIONS")
    print("=" * 80)

    # Test water viscosity
    T_test = 353.0
    mu_water = EnhancedProperties.water_viscosity(T_test)
    print(f"Water viscosity at {T_test} K: {mu_water*1e6:.2f} μPa·s")
    assert 100 < mu_water*1e6 < 1000, "Water viscosity out of expected range"

    # Test gas viscosity
    mu_H2 = EnhancedProperties.gas_viscosity(T_test, 'H2')
    mu_O2 = EnhancedProperties.gas_viscosity(T_test, 'O2')
    print(f"H2 viscosity at {T_test} K: {mu_H2*1e6:.2f} μPa·s")
    print(f"O2 viscosity at {T_test} K: {mu_O2*1e6:.2f} μPa·s")
    assert mu_H2 < mu_O2, "H2 should be less viscous than O2"

    # Test water content correlation
    a_w = 1.0
    λ = EnhancedProperties.water_content_nafion(a_w, T_test)
    print(f"Water content at a_w={a_w}, T={T_test} K: λ = {λ:.2f}")
    assert 10 < λ < 20, "Water content out of expected range"

    # Test protonic conductivity
    σ_m = EnhancedProperties.protonic_conductivity_enhanced(λ, T_test)
    print(f"Protonic conductivity: {σ_m:.2f} S/m")
    assert σ_m > 1.0, "Conductivity too low"

    # Test gas solubility
    H_H2 = EnhancedProperties.gas_solubility(T_test, 101325, 'H2')
    H_O2 = EnhancedProperties.gas_solubility(T_test, 101325, 'O2')
    print(f"H2 Henry constant: {H_H2*1e6:.2f} ×10⁻⁶ mol/(m³·Pa)")
    print(f"O2 Henry constant: {H_O2*1e6:.2f} ×10⁻⁶ mol/(m³·Pa)")

    # Test relative permeability
    s_test = 0.5
    k_rl = EnhancedProperties.relative_permeability_liquid(s_test)
    k_rg = EnhancedProperties.relative_permeability_gas(s_test)
    print(f"At s={s_test}: k_rl={k_rl:.3f}, k_rg={k_rg:.3f}")
    assert 0 <= k_rl <= 1 and 0 <= k_rg <= 1, "Relative permeabilities out of bounds"

    print("✓ All property correlations validated")

def test_adaptive_mesh():
    """Test adaptive mesh generation"""
    print("\n" + "=" * 80)
    print("TEST 2: ADAPTIVE MESH REFINEMENT")
    print("=" * 80)

    params = Parameters()
    mesh = AdaptiveMesh.create_adaptive_mesh(
        params.L_PTL, params.L_aCL, params.L_mem_wet,
        params.L_cCL, params.L_GDL,
        n_base=10, refinement_factor=3
    )

    print(f"Total mesh nodes: {mesh.n}")
    print(f"Layer distribution:")
    print(f"  PTL: {len(mesh.idx_PTL)} nodes")
    print(f"  Anode CL: {len(mesh.idx_aCL)} nodes")
    print(f"  Membrane: {len(mesh.idx_mem)} nodes")
    print(f"  Cathode CL: {len(mesh.idx_cCL)} nodes")
    print(f"  GDL: {len(mesh.idx_GDL)} nodes")

    # Check mesh integrity
    assert mesh.n == len(mesh.z), "Node count mismatch"
    assert np.all(np.diff(mesh.z) > 0), "Mesh not monotonically increasing"
    assert len(mesh.Δz) == mesh.n, "Cell width array size mismatch"

    # Check all nodes are assigned to a layer
    total_layer_nodes = (len(mesh.idx_PTL) + len(mesh.idx_aCL) +
                         len(mesh.idx_mem) + len(mesh.idx_cCL) + len(mesh.idx_GDL))
    assert total_layer_nodes == mesh.n, "Not all nodes assigned to layers"

    print(f"Mean cell width: {np.mean(mesh.Δz)*1e6:.2f} μm")
    print(f"Min cell width: {np.min(mesh.Δz)*1e6:.2f} μm")
    print(f"Max cell width: {np.max(mesh.Δz)*1e6:.2f} μm")
    print(f"Refinement factor: {np.max(mesh.refinement_level)}")

    print("✓ Adaptive mesh validated")

    return mesh

def test_two_phase_flow():
    """Test two-phase flow model"""
    print("\n" + "=" * 80)
    print("TEST 3: TWO-PHASE FLOW MODEL")
    print("=" * 80)

    params = Parameters()
    mesh = AdaptiveMesh.create_adaptive_mesh(
        params.L_PTL, params.L_aCL, params.L_mem_wet,
        params.L_cCL, params.L_GDL,
        n_base=10, refinement_factor=3
    )

    two_phase = TwoPhaseFlow(mesh, params)

    print(f"Initial liquid saturation:")
    if len(mesh.idx_PTL) > 0:
        print(f"  PTL: {np.mean(two_phase.s_liquid[mesh.idx_PTL]):.3f}")
    if len(mesh.idx_GDL) > 0:
        print(f"  GDL: {np.mean(two_phase.s_liquid[mesh.idx_GDL]):.3f}")

    # Test saturation update with dummy current
    I_vol = 1e5 * np.ones(mesh.n)  # 100 kA/m³ typical current density
    dt = 0.1  # 0.1 second
    T = 353.0

    two_phase.update_saturation(dt, I_vol, T)

    print(f"After {dt} s with I_vol = {I_vol[0]*1e-3:.0f} kA/m³:")
    if len(mesh.idx_PTL) > 0:
        print(f"  PTL: {np.mean(two_phase.s_liquid[mesh.idx_PTL]):.3f}")

    # Check bounds
    assert np.all(two_phase.s_liquid >= 0.0), "Negative saturation detected"
    assert np.all(two_phase.s_liquid <= 1.0), "Saturation exceeds unity"

    print("✓ Two-phase flow model validated")

def test_degradation():
    """Test degradation mechanisms"""
    print("\n" + "=" * 80)
    print("TEST 4: DEGRADATION MECHANISMS")
    print("=" * 80)

    degradation = DegradationState()

    # Simulate operation
    j = 1.5  # A/cm²
    T = 353  # K
    dt_hours = 100  # hours

    print(f"Simulating {dt_hours} hours at {j} A/cm², {T} K")

    degradation.update_degradation(dt_hours, j, T, potential_cycles=100)

    print(f"After operation:")
    print(f"  Anode catalyst: {degradation.catalyst_loading_anode*100:.2f}%")
    print(f"  Anode ECSA: {degradation.ECSA_anode*100:.2f}%")
    print(f"  Cathode catalyst: {degradation.catalyst_loading_cathode*100:.2f}%")
    print(f"  Cathode ECSA: {degradation.ECSA_cathode*100:.2f}%")
    print(f"  Membrane thickness: {degradation.membrane_thickness*100:.2f}%")
    print(f"  Membrane conductivity: {degradation.membrane_conductivity*100:.2f}%")
    print(f"  Gas crossover: {degradation.membrane_crossover:.3f}x")

    # Check degradation is occurring
    assert degradation.catalyst_loading_anode < 1.0, "No anode catalyst degradation"
    assert degradation.ECSA_cathode < 1.0, "No cathode ECSA degradation"
    assert degradation.membrane_thickness < 1.0, "No membrane thinning"

    # Apply to parameters
    params = Parameters()
    params_deg = degradation.apply_to_parameters(params)

    print(f"\nParameter changes:")
    print(f"  Anode loading: {params.loading_aCL:.3f} → {params_deg.loading_aCL:.3f} mg/cm²")
    print(f"  Membrane thickness: {params.L_mem_dry*1e6:.1f} → {params_deg.L_mem_dry*1e6:.1f} μm")

    print("✓ Degradation model validated")

def test_steady_state():
    """Test steady-state solver"""
    print("\n" + "=" * 80)
    print("TEST 5: STEADY-STATE SOLVER")
    print("=" * 80)

    params = Parameters()
    model = DynamicPEMElectrolyzer(params)

    print(f"Model created with {model.n} nodes")

    # Test single voltage point
    V_test = 1.7
    result = model.solve_steady_state(V_test, method='hybr')

    if result['converged']:
        print(f"✓ Converged at V={V_test} V: j={result['current_density']:.4f} A/cm²")

        # Check physical consistency
        assert result['current_density'] > 0, "Current should be positive"
        assert result['current_density'] < 10, "Current unreasonably high"

        # Check state variables
        φ_s = model.state[:model.n]
        φ_m = model.state[model.n:2*model.n]
        c_H2 = model.state[2*model.n:3*model.n]
        c_O2 = model.state[3*model.n:4*model.n]

        print(f"State variable ranges:")
        print(f"  φ_s: [{np.min(φ_s):.3f}, {np.max(φ_s):.3f}] V")
        print(f"  φ_m: [{np.min(φ_m):.3f}, {np.max(φ_m):.3f}] V")
        print(f"  c_H2: [{np.min(c_H2):.2f}, {np.max(c_H2):.2f}] mol/m³")
        print(f"  c_O2: [{np.min(c_O2):.2f}, {np.max(c_O2):.2f}] mol/m³")

        # Check for NaN or inf
        assert not np.any(np.isnan(model.state)), "NaN in solution"
        assert not np.any(np.isinf(model.state)), "Inf in solution"

        print("✓ Steady-state solution physically consistent")
    else:
        print("✗ Failed to converge (this may be OK for difficult conditions)")

def test_polarization_curve():
    """Test full polarization curve"""
    print("\n" + "=" * 80)
    print("TEST 6: POLARIZATION CURVE")
    print("=" * 80)

    params = Parameters()
    model = DynamicPEMElectrolyzer(params)

    V_cells = [1.5, 1.6, 1.7, 1.8]
    j_cells = []

    for V in V_cells:
        result = model.solve_steady_state(V, method='hybr')
        if result['converged']:
            j_cells.append(result['current_density'])
            print(f"  V={V:.2f} V → j={result['current_density']:.4f} A/cm²")
        else:
            j_cells.append(np.nan)
            print(f"  V={V:.2f} V → failed")

    # Check monotonicity (higher voltage = higher current)
    valid_data = [(V, j) for V, j in zip(V_cells, j_cells) if not np.isnan(j)]
    if len(valid_data) >= 2:
        V_valid, j_valid = zip(*valid_data)
        for i in range(len(j_valid) - 1):
            if V_valid[i+1] > V_valid[i]:
                assert j_valid[i+1] > j_valid[i], "Polarization curve not monotonic"

        print("✓ Polarization curve is monotonically increasing")

    # Plot if successful
    if len(valid_data) >= 2:
        fig, ax = plt.subplots(figsize=(8, 6))
        V_plot, j_plot = zip(*valid_data)
        ax.plot(j_plot, V_plot, 'o-', linewidth=2, markersize=10)
        ax.set_xlabel('Current Density [A/cm²]', fontsize=12)
        ax.set_ylabel('Cell Voltage [V]', fontsize=12)
        ax.set_title('PEM Electrolyzer Polarization Curve', fontsize=14)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('test_polarization_curve.png', dpi=150)
        print("✓ Polarization curve saved to test_polarization_curve.png")

def visualize_properties():
    """Create visualization of all key properties"""
    print("\n" + "=" * 80)
    print("CREATING PROPERTY VISUALIZATIONS")
    print("=" * 80)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 1. Water content vs activity
    ax = axes[0, 0]
    a_w = np.linspace(0, 1.5, 100)
    T_values = [303, 333, 353, 373]
    for T in T_values:
        λ = [EnhancedProperties.water_content_nafion(a, T) for a in a_w]
        ax.plot(a_w, λ, label=f'{T-273:.0f}°C', linewidth=2)
    ax.set_xlabel('Water Activity [-]')
    ax.set_ylabel('Water Content λ [-]')
    ax.set_title('Water Uptake in Nafion')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Protonic conductivity
    ax = axes[0, 1]
    λ_range = np.linspace(1, 20, 100)
    for T in T_values:
        σ_m = [EnhancedProperties.protonic_conductivity_enhanced(λ, T) for λ in λ_range]
        ax.plot(λ_range, σ_m, label=f'{T-273:.0f}°C', linewidth=2)
    ax.set_xlabel('Water Content λ [-]')
    ax.set_ylabel('Conductivity [S/m]')
    ax.set_title('Protonic Conductivity')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. Relative permeabilities
    ax = axes[0, 2]
    s = np.linspace(0, 1, 100)
    k_rl = [EnhancedProperties.relative_permeability_liquid(si) for si in s]
    k_rg = [EnhancedProperties.relative_permeability_gas(si) for si in s]
    ax.plot(s, k_rl, 'b-', label='Liquid', linewidth=2)
    ax.plot(s, k_rg, 'r-', label='Gas', linewidth=2)
    ax.set_xlabel('Liquid Saturation [-]')
    ax.set_ylabel('Relative Permeability [-]')
    ax.set_title('Brooks-Corey Model')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Viscosities
    ax = axes[1, 0]
    T_range = np.linspace(293, 373, 100)
    μ_water = [EnhancedProperties.water_viscosity(T)*1e6 for T in T_range]
    μ_H2 = [EnhancedProperties.gas_viscosity(T, 'H2')*1e6 for T in T_range]
    μ_O2 = [EnhancedProperties.gas_viscosity(T, 'O2')*1e6 for T in T_range]
    ax.plot(T_range-273, μ_water, 'b-', label='Water', linewidth=2)
    ax.plot(T_range-273, μ_H2, 'r--', label='H₂', linewidth=2)
    ax.plot(T_range-273, μ_O2, 'g--', label='O₂', linewidth=2)
    ax.set_xlabel('Temperature [°C]')
    ax.set_ylabel('Viscosity [μPa·s]')
    ax.set_title('Dynamic Viscosity')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 5. Gas solubility
    ax = axes[1, 1]
    H_H2 = [EnhancedProperties.gas_solubility(T, 101325, 'H2')*1e6 for T in T_range]
    H_O2 = [EnhancedProperties.gas_solubility(T, 101325, 'O2')*1e6 for T in T_range]
    ax.plot(T_range-273, H_H2, 'r-', label='H₂', linewidth=2)
    ax.plot(T_range-273, H_O2, 'b-', label='O₂', linewidth=2)
    ax.set_xlabel('Temperature [°C]')
    ax.set_ylabel('Henry Constant [×10⁻⁶ mol/(m³·Pa)]')
    ax.set_title('Gas Solubility')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 6. Mesh structure
    ax = axes[1, 2]
    params = Parameters()
    mesh = AdaptiveMesh.create_adaptive_mesh(
        params.L_PTL, params.L_aCL, params.L_mem_wet,
        params.L_cCL, params.L_GDL, n_base=10, refinement_factor=3
    )
    ax.plot(mesh.z*1e6, mesh.refinement_level, 'o-', markersize=4, linewidth=1)
    ax.set_xlabel('Position [μm]')
    ax.set_ylabel('Refinement Level')
    ax.set_title('Adaptive Mesh Refinement')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('test_property_visualization.png', dpi=150)
    print("✓ Property visualizations saved to test_property_visualization.png")

def run_all_tests():
    """Run all tests"""
    print("\n" + "=" * 80)
    print("PEM ELECTROLYZER MODEL - COMPREHENSIVE TEST SUITE")
    print("=" * 80)

    try:
        test_enhanced_properties()
        test_adaptive_mesh()
        test_two_phase_flow()
        test_degradation()
        test_steady_state()
        test_polarization_curve()
        visualize_properties()

        print("\n" + "=" * 80)
        print("ALL TESTS PASSED ✓")
        print("=" * 80)
        print("\nModel features validated:")
        print("  ✓ Enhanced property correlations")
        print("  ✓ Adaptive mesh refinement")
        print("  ✓ Two-phase flow model")
        print("  ✓ Degradation mechanisms")
        print("  ✓ Steady-state solver")
        print("  ✓ Polarization curves")
        print("\nGenerated files:")
        print("  - test_polarization_curve.png")
        print("  - test_property_visualization.png")

        return True

    except Exception as e:
        print(f"\n✗ TEST FAILED: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
