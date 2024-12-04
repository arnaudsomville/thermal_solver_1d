"""Definition of all the inputs of the model."""

from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np

@dataclass
class BoundaryConditions:
    """Define the boundary conditions."""
    boundary_temperature: float # T(x=e, t)
    alpha_p: float # Emissivity coefficient of the boundary material, value between 0 and 1
    incident_flux: float # Incident flux in W/m2
    conv_flux: float # Conv flux in W/m2

boltzmann_constant = 5.67e-8

@dataclass
class PhysicalParameters:
    """Physical parameters."""
    material_name: str
    ablation_temperature: float
    rho: float # Material density
    cp: float # Specific heat capacity
    k: float # Thermal conductivity

class Cell:
    """Define a single cell."""
    def __init__(self, initial_temperature: float, physical_parameters: PhysicalParameters)->None:
        """_summary_

        Args:
            initial_temperature (float): T0
            initial_temperature (PhysicalParameters): Physical parameters of the cell
        """
        self.temperature: dict[int, float] = {0 : initial_temperature}
        self.physical_parameters: PhysicalParameters = physical_parameters
    
    def update_temperature(self, new_temperature: float, time_id: int)->None:
        """Update the temperature at the current time step."""
        self.temperature[time_id] = new_temperature
    
    def print_temperature(self, delta_t: float):
        """Print the temperature."""
        time_steps = np.array(list(self.temperature.keys()))
        temperatures = list(self.temperature.values())
        plt.plot(time_steps*delta_t, temperatures, label=self.physical_parameters.material_name)
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (K)')
        plt.legend()
        plt.show()