import math

import pytest

from ka_from_tmax_kel import estimate_ka, tmax_from_rates


@pytest.mark.parametrize("ka, kel", [(1.2, 0.2), (0.1, 0.5), (0.5, 0.5), (0.500000001, 0.5)])
def test_round_trip(ka, kel):
    tmax = tmax_from_rates(ka, kel)
    result = estimate_ka(tmax, kel)
    assert math.isclose(result.ka, ka, rel_tol=2e-9)
    assert abs(result.residual) <= 2e-10 * tmax


def test_known_value():
    result = estimate_ka(tmax=1.5, kel=0.5)
    assert result.ka == pytest.approx(0.8668004944, rel=1e-9)


@pytest.mark.parametrize("tmax, kel", [(0, 0.5), (-1, 0.5), (1, 0), (float("nan"), 0.5)])
def test_rejects_invalid_inputs(tmax, kel):
    with pytest.raises(ValueError):
        estimate_ka(tmax, kel)
