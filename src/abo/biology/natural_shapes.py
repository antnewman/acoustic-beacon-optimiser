"""Published reflector dimensions from the literature.

Data sources:
- Marcgravia evenia: Simon et al. (2011), Science 333: 631-633
- Mucuna holtonii: von Helversen and von Helversen (1999)
- Simon et al. (2020), PNAS 117: 1367-1374
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NaturalReflector:
    """Published dimensions of a natural reflector.

    Parameters
    ----------
    species : str
        Plant species name.
    structure : str
        Morphological structure (e.g. 'dish-leaf', 'vexillum').
    diameter_mm : float
        Approximate diameter in mm.
    depth_mm : float
        Approximate depth in mm.
    source : str
        Literature citation.
    """

    species: str
    structure: str
    diameter_mm: float
    depth_mm: float
    source: str


MARCGRAVIA_EVENIA = NaturalReflector(
    species="Marcgravia evenia",
    structure="dish-leaf",
    diameter_mm=35.0,
    depth_mm=25.0,
    source="Simon et al. (2011), Science 333: 631-633",
)

MUCUNA_HOLTONII = NaturalReflector(
    species="Mucuna holtonii",
    structure="vexillum",
    diameter_mm=25.0,
    depth_mm=15.0,
    source="von Helversen and von Helversen (1999)",
)

# Simon et al. (2020) spherical cap configurations.
# Note: r=35 with d=49 is geometrically impossible (d > r) and excluded.
SIMON_2020_CONFIGS: list[dict[str, float]] = [
    {"radius_mm": r, "depth_mm": d}
    for r in [35.0, 50.0, 70.0]
    for d in [25.0, 30.0, 49.0]
    if d <= r
]
