
from td_rentee_atmo.mesher import Mesher, ThermalProtectionLayer
from td_rentee_atmo.solver import ThermalSolver
from td_rentee_atmo.solver_inputs import BoundaryConditions, PhysicalParameters


def test_global_solve():
    cork_params =PhysicalParameters(
        material_name='cork',
        ablation_temperature=1000,
        rho=120,
        cp = 1900,
        k = 0.04
    )

    silica_params =PhysicalParameters( 
        material_name='LI-900 tiles',
        ablation_temperature=10000,
        rho=144,
        cp = 1250,
        k = 0.02
    )

    layers = [
        ThermalProtectionLayer(width_mm=50, t0=20, physical_parameter_m=cork_params),
        ThermalProtectionLayer(width_mm=32, t0=30, physical_parameter_m=silica_params),
    ]

    delta_x = 0.5
    delta_t = 0.05
    mesher = Mesher()
    mesh = mesher.mesh(layers, delta_x)
    mesher.print_meshing()

    hypersonic_conditions = BoundaryConditions(
        boundary_temperature=800,
        alpha_p=0.8,
        incident_flux=30e3,
        conv_flux=5e3
    )

    solver = ThermalSolver(mesh, delta_x, delta_t, hypersonic_conditions)
    
    solver.solve(100)
    solver.print_solution()
