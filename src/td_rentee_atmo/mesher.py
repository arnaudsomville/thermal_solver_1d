"""Definition of a 1 dimensional mesher."""

from dataclasses import dataclass
import numpy as np
from td_rentee_atmo.solver_inputs import Cell, PhysicalParameters


@dataclass
class ThermalProtectionLayer:
    """Define the thermal protection."""
    width_mm: float
    t0: float
    physical_parameter_m: PhysicalParameters

class Mesher:
    """Meshes a line."""

    def __init__(self)->None:
        """Initialize empty Mesh."""
        self.meshing: dict[int, Cell] = {}
        self.delta_x: float = 0
        pass

    def set_meshing(self, meshing: dict[int, Cell], delta_x: float)->None:
        """_summary_

        Args:
            meshing (dict[int, Cell]): list of cells in mesh.
            delta_x (float): distance between cells in mm.
        """
        self.meshing = meshing
        self.delta_x = delta_x

    
    def mesh(self, layers: list[ThermalProtectionLayer], delta_x: float)->dict[int, Cell]:
        """Divides the thermal protection into a list of cells.
        
        Args:
            layers (list[ThermalProtectionLayer]): list containing layers data.
            delta_x (float): distance between two cells in mm.
        
        Returns:
            dict[int, Cell]: list of cells.
        """

        @dataclass
        class LayerBoundary:
            layer_start: float
            layer_end: float

        thermal_protection_lenght: float = 0
        last_boundary_end: float = 0
        layer_boundaries: dict[int, LayerBoundary] = {}

        for i, layer in enumerate(layers):
            layer_boundaries[i] = LayerBoundary(
                layer_start = last_boundary_end,
                layer_end = layer.width_mm
            )
            thermal_protection_lenght += layer.width_mm
            last_boundary_end = layer.width_mm

        nb_cells = int(np.ceil(thermal_protection_lenght / delta_x) + 1)
        
        cell_dict: dict[int, Cell] = {}
        for i in range(nb_cells):
            for bound_id, bounds in layer_boundaries.items():
                if i*delta_x >= bounds.layer_start and i*delta_x < bounds.layer_end:
                    break
            cell_dict[i] = Cell(
                initial_temperature=layers[bound_id].t0,
                physical_parameters=layers[bound_id].physical_parameter_m,
            )
        self.meshing = cell_dict
        self.delta_x = delta_x

        return cell_dict
    
    def print_meshing(self)->None:
        """Print Meshing infos."""
        for i, cell in self.meshing.items():
            print(f"From {i*self.delta_x}mm : {cell.physical_parameters.material_name}")
    
if __name__ == '__main__':
    cork_params =PhysicalParameters(
        material_name='cork',
        ablation_temperature=10000,
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
        ThermalProtectionLayer(width_mm=32, t0=30, physical_parameter_m=silica_params),
        ThermalProtectionLayer(width_mm=50, t0=20, physical_parameter_m=cork_params),
    ]

    delta_x = 10
    mesher = Mesher()
    mesh = mesher.mesh(layers, delta_x)
    mesher.print_meshing()
