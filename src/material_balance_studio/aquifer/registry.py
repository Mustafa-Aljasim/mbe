"""Aquifer model selections create new immutable configurations."""
from .none import NoAquifer
from .pot import PotAquifer
from .schilthuis import SchilthuisAquifer
from .fetkovich import FetkovichAquifer
from .carter_tracy import CarterTracyAquifer
from .veh import VanEverdingenHurstAquifer
from .modified_veh import ModifiedVanEverdingenHurstAquifer

MODELS = {"None": NoAquifer, "Pot": PotAquifer, "Schilthuis": SchilthuisAquifer,
          "Fetkovich": FetkovichAquifer, "Carter-Tracy": CarterTracyAquifer,
          "Van Everdingen-Hurst": VanEverdingenHurstAquifer,
          "Modified Van Everdingen-Hurst": ModifiedVanEverdingenHurstAquifer}
