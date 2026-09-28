# Ka from Tmax and Kel

A small, dependency-free Python package for estimating the first-order
absorption rate constant **Ka** from a reported **Tmax** and elimination rate
constant **Kel**.

**GitHub Topics:** `pharmacokinetics`, `absorption-rate-constant`, `tmax`,
`drug-absorption`

The calculation uses the standard one-compartment model with first-order oral
absorption and elimination:

```text
Tmax = ln(Ka / Kel) / (Ka - Kel)
```

The Python implementation replaces the fixed-step loop in the original MATLAB
version with a bracketed bisection solver. This gives a reproducible result with
a quantified residual and also handles the limiting case `Ka = Kel`.

## Assumptions

- One-compartment disposition
- First-order absorption and elimination
- No absorption lag time
- `Ka` and `Kel` are positive
- `Tmax` and the rate constants use compatible units (for example, hours and
  `h⁻¹`)

This is a modelling/research utility, not medical advice. A literature Tmax may
be affected by sampling design, formulation, food, population variability, and
model misspecification; report those sources and assumptions with the estimate.

## Installation

Requires Python 3.9 or newer.

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e '.[test]'
pytest
```

To enable diagnostic plots:

```bash
python -m pip install -e '.[plot]'
```

## Command line

```bash
ka-from-tmax-kel 1.5 0.5
```

Example output:

```text
Ka = 0.86680049
Residual Tmax = ...
Iterations = ...
```

Save a plot of the solution:

```bash
ka-from-tmax-kel 1.5 0.5 --plot ka_estimate.png
```

## Python API

```python
from ka_from_tmax_kel import estimate_ka

result = estimate_ka(tmax=1.5, kel=0.5)
print(result.ka)        # 0.866800...
print(result.residual) # calculated Tmax - supplied Tmax
```

## Repository contents

- `src/ka_from_tmax_kel/`: maintained Python package
- `tests/`: automated tests

## Scientific note

For fixed positive `Kel`, the equation is monotonic in `Ka`, so there is one
positive solution. If `Tmax < 1/Kel`, then `Ka > Kel`; if `Tmax > 1/Kel`, then
`Ka < Kel` (a flip-flop-like solution). At `Tmax = 1/Kel`, the limiting solution
is `Ka = Kel`.

## License

MIT. See [LICENSE](LICENSE).
