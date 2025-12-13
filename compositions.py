# compositions.py
"""
Material property database for meteoroid compositions.

This module provides thermophysical properties for common meteoroid types:
- Iron: Fe-Ni metal (iron meteorites, ~5% of falls)
- Stone: Silicate rocks (ordinary chondrites, ~85% of falls)
- Ice: Water ice (cometary material, volatile)
- Carbonaceous: Carbon-rich material (carbonaceous chondrites, ~4% of falls)

Properties include:
- Thermal: density, specific heat, thermal conductivity, emissivity
- Phase change: fusion/vaporization temperatures, latent heats
- Mechanical: tensile strength, drag coefficient

References:
- Ceplecha et al. (1998) "Meteor Phenomena and Bodies"
- Popova et al. (2011) "Very Low Strengths of Interplanetary Meteoroids"
- Standard thermophysical property databases (NIST, etc.)
"""


class Compositions:
    def __init__(self):
        self.library = {
            # Iron meteorites (Fe-Ni alloy, ~7-8% Ni typical)
            # High density, high thermal conductivity, resistant to ablation
            "iron": {
                "density": 7800.0,           # kg/m³ (Fe-Ni alloy)
                "cp": 450.0,                 # J/(kg·K) (specific heat)
                "thermal_k": 80.0,           # W/(m·K) (thermal conductivity)
                "emissivity": 0.75,          # Oxidized iron surface
                "fusion_T": 1809.0,          # K (melting point of iron)
                "vapor_T": 3134.0,           # K (boiling point at 1 atm)
                "L_fusion": 2.7e5,           # J/kg (latent heat of fusion)
                "L_vapor": 6.3e6,            # J/kg (latent heat of vaporization)
                "tensile": 1.0e8,            # Pa (tensile strength)
                "Cd": 1.0,                   # Drag coefficient (sphere)
                "volatile": False,
                "description": "Iron-nickel meteorite (siderite)"
            },

            # Stony meteorites (silicate minerals: olivine, pyroxene)
            # Most common type, moderate strength
            "stone": {
                "density": 3500.0,           # kg/m³ (typical ordinary chondrite)
                "cp": 800.0,                 # J/(kg·K)
                "thermal_k": 2.0,            # W/(m·K) (silicates, low conductivity)
                "emissivity": 0.90,          # High emissivity silicates
                "fusion_T": 1700.0,          # K (silicate melting, solidus)
                "vapor_T": 2800.0,           # K (silicate vaporization)
                "L_fusion": 5.0e5,           # J/kg
                "L_vapor": 5.0e6,            # J/kg
                "tensile": 2.0e7,            # Pa (fractured rock strength)
                "Cd": 1.0,
                "volatile": False,
                "description": "Ordinary chondrite (stony meteorite)"
            },

            # Water ice (cometary material)
            # Low density, sublimates at low pressures
            "ice": {
                "density": 917.0,            # kg/m³ (ice Ih at 0°C)
                "cp": 2090.0,                # J/(kg·K) (ice at 0°C)
                "thermal_k": 2.2,            # W/(m·K)
                "emissivity": 0.95,          # High emissivity
                "fusion_T": 273.15,          # K (melting point at 1 atm)
                "vapor_T": 373.15,           # K (boiling at 1 atm; note: sublimes at P < 611 Pa)
                "L_fusion": 3.34e5,          # J/kg (latent heat of fusion)
                "L_vapor": 2.26e6,           # J/kg (vaporization) / 2.83e6 (sublimation)
                "L_sublimation": 2.83e6,     # J/kg (direct solid→vapor)
                "tensile": 1.0e6,            # Pa (ice fracture strength)
                "Cd": 0.9,                   # Slightly streamlined due to ablation shape
                "volatile": True,
                "triple_point_P": 611.657,   # Pa (triple point pressure)
                "description": "Water ice (cometary nucleus material)"
            },

            # Carbonaceous chondrites (carbon-rich, porous, primitive)
            # Contains organics, hydrated minerals, volatiles
            # Low strength, high porosity (10-40%)
            "carbonaceous": {
                "density": 2200.0,           # kg/m³ (lower due to porosity)
                "cp": 900.0,                 # J/(kg·K)
                "thermal_k": 1.0,            # W/(m·K) (porous, low conductivity)
                "emissivity": 0.90,          # Dark, carbonaceous surface
                "fusion_T": 1500.0,          # K (lower due to volatile content)
                "vapor_T": 2200.0,           # K (organics volatilize earlier)
                "L_fusion": 4.0e5,           # J/kg
                "L_vapor": 4.0e6,            # J/kg
                "tensile": 1.0e7,            # Pa (weak, friable)
                "Cd": 1.1,                   # Higher drag due to irregular shape
                "volatile": True,
                "porosity": 0.25,            # Typical porosity fraction
                "description": "Carbonaceous chondrite (primitive, volatile-rich)"
            },
        }

    def get(self, name: str) -> dict:
        """
        Retrieve material properties.

        Parameters
        ----------
        name : str
            Material name (case-insensitive)

        Returns
        -------
        dict
            Copy of material properties

        Raises
        ------
        ValueError
            If material not in library
        """
        key = name.lower()
        if key not in self.library:
            available = ", ".join(self.library.keys())
            raise ValueError(f"Unknown composition '{name}'. Available: {available}")
        return self.library[key].copy()

    def list_materials(self) -> list:
        """Return list of available material names."""
        return list(self.library.keys())

    def get_property(self, name: str, prop: str):
        """
        Get a specific property for a material.

        Parameters
        ----------
        name : str
            Material name
        prop : str
            Property key (e.g., 'density', 'fusion_T')

        Returns
        -------
        Property value or None if not found
        """
        mat = self.get(name)
        return mat.get(prop)
