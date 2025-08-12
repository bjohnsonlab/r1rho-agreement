## Eric's MATLAB Data

Load MATLAB and run `eric.m`

### Parameters used

* νC = 99.8 MHz => B₀ ≈ 396.8 MHz
* νR = 10 kHz
* δ = -37 ppm => Δσ = -55.5 ppm
* S² = 0.3295
* τc:
  - 37°C: τc = 56.6397 ns => log₁₀(τc) = -7.246879
  - 57°C: τc = 350.217 μs => log₁₀(τc) = -3.455660
  - 67°C: τc = 147.259 μs => log₁₀(τc) = -3.831919
  - 77°C: τc = 104.318 μs => log₁₀(τc) = -3.981640

## RING Data

* In `File > Preferences...`:
  - Set `Limits > Reference Field` to 396.8
  - Set `ModelFree > C CSA` to 55.5
* Open the `Simulate` tab
* Select `SSR1Rho` and `CSA`
* Set the `S²` slider to 0.3295
* Set the `νR` slider to 10.0
* For each temparature (`T`) and `log₁₀(τc)` pairing above:
  - Set the `log₁₀(τc)` slider
  - Click `File > Save XY Chart Data`.
  - Save the text file (`data/ring_{T}.txt`)
