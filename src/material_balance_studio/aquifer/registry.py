"""Only Phase 3A models. Selections create new immutable configurations."""
from .none import NoAquifer
from .pot import PotAquifer
from .schilthuis import SchilthuisAquifer
from .fetkovich import FetkovichAquifer

MODELS = {"None": NoAquifer, "Pot": PotAquifer, "Schilthuis": SchilthuisAquifer, "Fetkovich": FetkovichAquifer}
