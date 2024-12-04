"""File containong the definition of the 1D solver."""


import numpy as np
from td_rentee_atmo.mesher import Mesher, ThermalProtectionLayer
from td_rentee_atmo.solver_inputs import BoundaryConditions, Cell, PhysicalParameters, boltzmann_constant


class ThermalSolver:
    """1D thermal solver."""
    def __init__(self,meshing: dict[int, Cell], delta_x: float, delta_t: float, boundary_conditions: BoundaryConditions)->None:
        """Initialize solver.
        
        Args:
            delta_x (float): x step.
            delta_t (float): y step.
        """
        self.input_meshing = meshing
        self.delta_x = delta_x
        self.delta_t = delta_t
        self.cell_stability_coefficients = self.compute_s()
        self.verify_cell_stability_coefficients()
        
        self.boundary_conditions = boundary_conditions

    def compute_s(self)->dict[int, float]:
        """Compute the dictionnary associating to each cell of the meshing the stability value. (should be < 1/2)  

        Returns:
            dict[int, float]: Dictionnary associating to each cell of the meshing the stability.
        """
        cell_stability_coefficients = {}
        for cell_id, cell in self.input_meshing.items():
            cell_stability_coefficients[cell_id] =cell.physical_parameters.k * self.delta_t / ((self.delta_x*1e-3) ** 2 * cell.physical_parameters.rho* cell.physical_parameters.cp)
        return cell_stability_coefficients
    
    def verify_cell_stability_coefficients(self)->None:
        """Verify that all cell stability coefficients are < 1/2."""
        for cell_id, s in self.cell_stability_coefficients.items():
            if s >= 0.5:
                raise ValueError(f'Error: Cell {cell_id} has a stability coefficient {s} which is greater than 1/2.')
        print("Meshing stability verified successfully !")

    def solve(self, duration_s: float)->dict[int, Cell]:
        """Solve the Thermal problem.

        Args:
            duration_s (float): Amount of time of the simulation.

        Returns:
            dict[int, Cell]: Dictionnary containing the meshing for each instant.

        """
        self.temporal_meshing_solution = self.input_meshing.copy() # Reinitialize Meshing Solution
        first_non_ablated_cell = 0
        for t in range(first_non_ablated_cell+1, int(np.ceil(duration_s / self.delta_t)) + 1):
            precedent_meshing_temp = {cell_id : cell.temperature[t-1] for cell_id, cell in self.input_meshing.items()}
            new_meshing_temp = precedent_meshing_temp.copy()
            # Dirichlet condition
            flux_conv = self.boundary_conditions.conv_flux
            incident_flux = self.boundary_conditions.alpha_p*self.boundary_conditions.incident_flux
            flux_rad = self.boundary_conditions.alpha_p*boltzmann_constant*new_meshing_temp[first_non_ablated_cell]**4
            if flux_rad >= incident_flux:
                flux_rad = incident_flux
            new_meshing_temp[first_non_ablated_cell] = precedent_meshing_temp[first_non_ablated_cell+1] - 2*self.delta_x*1e-3*(flux_rad - incident_flux- flux_conv)/self.input_meshing[first_non_ablated_cell].physical_parameters.k
            
            # Neumann condition
            new_meshing_temp[len(new_meshing_temp)-1] = self.boundary_conditions.boundary_temperature
            
            for cell_id in range(first_non_ablated_cell+1, len(new_meshing_temp)-1):
                new_meshing_temp[cell_id] = self.cell_stability_coefficients[cell_id]*(precedent_meshing_temp[cell_id-1] -2*precedent_meshing_temp[cell_id]+ precedent_meshing_temp[cell_id+1]) + precedent_meshing_temp[cell_id]

            # Ablation in the side opposite of the boundary
            for cell_id in range(first_non_ablated_cell, len(new_meshing_temp)-1):   
                if new_meshing_temp[cell_id] >= self.temporal_meshing_solution[cell_id].physical_parameters.ablation_temperature:
                   new_meshing_temp[cell_id] = -1
                   first_non_ablated_cell += 1
                else:
                    break


            if first_non_ablated_cell > 0:
                for cell_id in range(0, first_non_ablated_cell):
                    new_meshing_temp[cell_id] = -1

            for cell_id, cell in self.temporal_meshing_solution.items():
                cell.temperature[t] = new_meshing_temp[cell_id]
            
            if first_non_ablated_cell >= len(new_meshing_temp)-2:
                print(f"Thermal protection completely burned before the complete simulation. It last {t*self.delta_t}s")
                return self.temporal_meshing_solution
            

        return self.temporal_meshing_solution

    def print_solution(self)->None:
        """Print the meshing's solution."""
        for t in self.temporal_meshing_solution[0].temperature:
            print(f"Time {t*self.delta_t}:")
            print(''.join([f'|{cell_data.temperature[t]}°C' for cell_data in self.temporal_meshing_solution.values()]))
            print("------------------------------------------------------------------------")

if __name__ == '__main__': # pragma: no cover
    cork_params =PhysicalParameters(
        material_name='cork',
        ablation_temperature=485,
        rho=120,
        cp = 1900,
        k = 0.04
    )

    silica_params =PhysicalParameters( 
        material_name='LI-900 tiles',
        ablation_temperature=1500,
        rho=144,
        cp = 1250,
        k = 0.02
    )

    layers = [
        ThermalProtectionLayer(width_mm=32, t0=30, physical_parameter_m=silica_params),
        ThermalProtectionLayer(width_mm=50, t0=20, physical_parameter_m=cork_params),
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
    #solver.temporal_meshing_solution[int(len(solver.temporal_meshing_solution)/2)].print_temperature(delta_t)
    solver.temporal_meshing_solution[0].print_temperature(delta_t)