"""3D mesh rendering using PyVista.

Renders reflector surface meshes with optional field overlays
(e.g. surface pressure magnitude).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    pass


def plot_mesh(
    grid: object,
    scalars: NDArray[np.float64] | None = None,
    scalar_label: str = "Pressure (Pa)",
) -> None:
    """Render a reflector mesh in 3D.

    Parameters
    ----------
    grid : bempp.api.Grid
        Surface mesh.
    scalars : NDArray[np.float64], optional
        Scalar field on mesh vertices for colouring.
    scalar_label : str
        Label for the colour bar.
    """
    raise NotImplementedError
