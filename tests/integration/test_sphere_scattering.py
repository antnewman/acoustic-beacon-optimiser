"""Integration test: validate BEM against analytical Mie series for a rigid sphere.

At He = 5, 10, 20, 30, the BEM-computed target strength should agree
with the Mie series solution to within 0.5 dB.
"""

from __future__ import annotations

import pytest


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.skip(reason="BEM solver not yet implemented")
class TestMieSeriesValidation:
    def test_he_5(self) -> None:
        pass

    def test_he_10(self) -> None:
        pass

    def test_he_20(self) -> None:
        pass

    def test_he_30(self) -> None:
        pass
