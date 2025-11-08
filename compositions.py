class Compositions:
    def __init__(self):
        self.library = {
            "iron": {
                "density": 7800.0, "cp": 450.0, "emissivity": 0.75,
                "fusion_T": 1809.0, "vapor_T": 3134.0,
                "L_fusion": 2.7e5, "L_vapor": 6.3e6,
                "tensile": 1.0e8, "Cd": 1.0,
                "volatile": False,
            },
            "stone": {
                "density": 3500.0, "cp": 800.0, "emissivity": 0.90,
                "fusion_T": 1700.0, "vapor_T": 2500.0,
                "L_fusion": 5.0e5, "L_vapor": 5.0e6,
                "tensile": 2.0e7, "Cd": 1.0,
                "volatile": False,
            },
            "ice": {
                "density": 900.0, "cp": 2100.0, "emissivity": 0.95,
                "fusion_T": 273.15, "vapor_T": 373.15,   # 1 atm reference
                "L_fusion": 3.3e5, "L_vapor": 2.6e6,
                "tensile": 1.0e6, "Cd": 0.9,
                "volatile": True,
            },
            "carbonaceous": {
                "density": 2200.0, "cp": 900.0, "emissivity": 0.90,
                "fusion_T": 1500.0, "vapor_T": 2200.0,
                "L_fusion": 4.0e5, "L_vapor": 4.0e6,
                "tensile": 1.0e7, "Cd": 1.1,
                "volatile": True,
            },
        }

    def get(self, name: str):
        key = name.lower()
        if key not in self.library:
            raise ValueError(f"Unknown composition '{name}'")
        return self.library[key].copy()