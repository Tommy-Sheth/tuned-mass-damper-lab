# Tuned Mass Damper Lab

Interactive simulation for our PH 217 group project *Stabilizing a Tower: the Taipei 101 Tuned Mass Damper*
(Rudra, Vineet, Riya Shinde; Prof. Prashant Kumar).

It models the first sway mode of Taipei 101 as one mass on a spring, with a 660 t pendulum damper hanging inside it, and
shows how the damper reduces sway in wind and in earthquakes.

## Run it

It is one self-contained file with no libraries and no internet needed.

- Download or clone the repository and open `index.html` in any browser, or
- (live version) turn on GitHub Pages: *Settings -> Pages -> Deploy from a branch -> main -> / (root)*,
  then the page is at `https://<your-username>.github.io/<repo-name>/`.

Pick a load (wind storm, steady push, frequency sweep, single gust, release, El Centro, long-period shaking, near-fault pulse),
press **Play**, and change the sliders for the building and the damper. The *Tuning* part draws the map of average sway against
tuning and damping and finds the best settings for four different goals.

## Is it correct?

Two independent checks are included.

1. **The page checks itself.** On load it runs 43 comparisons against reference values from our Python model
   (fixed points, peaks, average sway, wind, time stepping, energy balance, earthquakes, optimisers). To see them in the page, open
   the "checks" section at the bottom. To run them headless:

   ```
   pip install playwright && playwright install chromium
   python verify/run_page_checks.py
   ```
   Expected result (saved in `verify/expected_page_checks.txt`): *"The page and the Python model agree on all 43 checks."*

2. **A separate Python model** (`verify/tmdcheck.py`, written from the equations only) reproduces the main numbers of the report:

   ```
   pip install numpy scipy
   python verify/python_check.py
   ```
   Expected output is in `verify/expected_python_output.txt`, for example Design B `f = 0.9868`, `zeta_d = 0.0641`, `L = 11.34 m`,
   average sway ratio `0.584`, and 7.15 -> 4.72 milli-g in 50 random storms (a 34% cut).

## Files

| File | What it is |
|---|---|
| `index.html` | the whole simulation (physics, plots, animation, controls, embedded earthquake records and fonts) |
| `verify/tmdcheck.py` | independent Python model |
| `verify/python_check.py` | prints the main numbers of the report |
| `verify/run_page_checks.py` | runs the page's own 43 checks in a headless browser |
| `verify/expected_*.txt` | the outputs we got when we ran the two scripts |

## AI use

The simulation page was first coded by Rudra; Claude (Anthropic) was used to structure it (physics separated from the page, and the
self-checks against Python). `verify/tmdcheck.py` and the check scripts were written with Claude from the equations in the report.
