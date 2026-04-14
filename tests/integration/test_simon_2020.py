"""Integration test: validate BEM against Simon et al. (2020) FEM results.

Compares BEM-computed impulse responses and target strength for spherical
cap reflectors with r = {35, 50, 70} mm and d = {25, 30, 49} mm against
published FEM results from Simon et al. (2020), PNAS 117: 1367-1374.
"""

from __future__ import annotations

import pytest


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.skip(reason="BEM solver not yet implemented")
class TestSimon2020Validation:
    def test_spectral_notch_positions(self) -> None:
        pass

    def test_target_strength_agreement(self) -> None:
        pass

    def test_impulse_response_peaks(self) -> None:
        pass
