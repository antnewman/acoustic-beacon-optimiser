"""Integration test: validate BEM against analytical Mie series for a rigid sphere.

At He = 5 and He = 10, the BEM-computed monostatic target strength should
agree with the Mie series solution to within 0.5 dB when using a
sufficiently fine mesh (at least 6 elements per wavelength).

Higher Helmholtz numbers (He = 20, 30) require very fine meshes and are
marked as slow tests.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

sys.modules["bempp"] = __import__("bempp_cl")
sys.modules["bempp.api"] = __import__("bempp_cl.api", fromlist=["api"])
import bempp.api  # noqa: E402

from abo.acoustics.bem_solver import (  # noqa: E402
    scattered_field,
    solve_scattering,
)
from tests.analytical.mie_series import mie_backscatter_ts  # noqa: E402

SPHERE_RADIUS = 1.0
SPEED_OF_SOUND = 343.0
FAR_FIELD_RANGE = 10.0
TOLERANCE_DB = 0.5


def _frequency_for_he(he: float) -> float:
    """Compute frequency that gives the specified Helmholtz number."""
    return he * SPEED_OF_SOUND / (2.0 * np.pi * SPHERE_RADIUS)


def _mesh_level_for_he(he: float) -> int:
    """Choose bempp regular_sphere refinement level for adequate resolution.

    Targets at least 6 elements per wavelength on the sphere surface.
    """
    if he <= 5:
        return 4  # 2048 elements
    if he <= 12:
        return 5  # 8192 elements
    if he <= 25:
        return 6  # 32768 elements
    return 7  # 131072 elements


def _bem_backscatter_ts(he: float) -> float:
    """Compute BEM monostatic TS for a rigid unit sphere at given He."""
    f = _frequency_for_he(he)
    level = _mesh_level_for_he(he)
    grid = bempp.api.shapes.regular_sphere(level)

    direction = np.array([0.0, 0.0, -1.0])
    total = solve_scattering(grid, f, direction, SPEED_OF_SOUND)

    backscatter_pt = np.array([[0.0], [0.0], [FAR_FIELD_RANGE]])
    p_scat = scattered_field(
        grid, total, backscatter_pt, f, direction, SPEED_OF_SOUND,
    )

    return float(
        10.0 * np.log10(np.abs(p_scat[0]) ** 2 * FAR_FIELD_RANGE**2),
    )


@pytest.mark.integration
class TestMieSeriesValidation:
    def test_he_5(self) -> None:
        he = 5
        f = _frequency_for_he(he)
        ts_mie = mie_backscatter_ts(
            SPHERE_RADIUS, f, SPEED_OF_SOUND, FAR_FIELD_RANGE,
        )
        ts_bem = _bem_backscatter_ts(he)
        assert abs(ts_bem - ts_mie) < TOLERANCE_DB, (
            f"He={he}: BEM={ts_bem:.2f} dB, Mie={ts_mie:.2f} dB, "
            f"diff={abs(ts_bem - ts_mie):.3f} dB"
        )

    def test_he_10(self) -> None:
        he = 10
        f = _frequency_for_he(he)
        ts_mie = mie_backscatter_ts(
            SPHERE_RADIUS, f, SPEED_OF_SOUND, FAR_FIELD_RANGE,
        )
        ts_bem = _bem_backscatter_ts(he)
        assert abs(ts_bem - ts_mie) < TOLERANCE_DB, (
            f"He={he}: BEM={ts_bem:.2f} dB, Mie={ts_mie:.2f} dB, "
            f"diff={abs(ts_bem - ts_mie):.3f} dB"
        )

    @pytest.mark.slow
    def test_he_20(self) -> None:
        he = 20
        f = _frequency_for_he(he)
        ts_mie = mie_backscatter_ts(
            SPHERE_RADIUS, f, SPEED_OF_SOUND, FAR_FIELD_RANGE,
        )
        ts_bem = _bem_backscatter_ts(he)
        assert abs(ts_bem - ts_mie) < TOLERANCE_DB, (
            f"He={he}: BEM={ts_bem:.2f} dB, Mie={ts_mie:.2f} dB, "
            f"diff={abs(ts_bem - ts_mie):.3f} dB"
        )

    @pytest.mark.slow
    def test_he_30(self) -> None:
        he = 30
        f = _frequency_for_he(he)
        ts_mie = mie_backscatter_ts(
            SPHERE_RADIUS, f, SPEED_OF_SOUND, FAR_FIELD_RANGE,
        )
        ts_bem = _bem_backscatter_ts(he)
        assert abs(ts_bem - ts_mie) < TOLERANCE_DB, (
            f"He={he}: BEM={ts_bem:.2f} dB, Mie={ts_mie:.2f} dB, "
            f"diff={abs(ts_bem - ts_mie):.3f} dB"
        )
