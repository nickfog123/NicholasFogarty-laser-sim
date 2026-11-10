# Sound-Driven Laser Pattern Simulator (MATLAB)

`LaserSpeakerSimulator.m` is a single, self-contained MATLAB app. It simulates the speaker, tape, spring and mirror apparatus in the photographs and shows **one continuously animated red laser trajectory** on a dark projection screen. Physics produces the trajectory in this order: tone, then amplifier, then both speakers, then the taped 1-inch spring, then mirror tilt, then exact ray reflection.

> All mechanical constants are **illustrative estimates** chosen to give recognisable patterns. The app is not calibrated, so do not treat it as an exact reproduction of this specific apparatus.

---

## 1. Running the program

1. Requires **MATLAB R2020b or newer** (`uifigure`, `uigridlayout`, `WordWrap`). No toolboxes are needed. `ode45` is part of core MATLAB.
2. Put `LaserSpeakerSimulator.m` in your working folder (or anywhere on the MATLAB path).
3. Start the GUI:
   ```matlab
   LaserSpeakerSimulator
   ```
   Press **Play**. The default is A4 = 440 Hz, a sine wave at 50 % volume, which gives a narrow tilted ellipse. Drag the frequency slider while the app is playing to see lines, ellipses and other figures.
4. Run the numerical checks without the GUI:
   ```matlab
   R = LaserSpeakerSimulator('verify');   % prints a PASS/FAIL table, returns a struct array
   ```
5. Closing the window stops and deletes the animation timer.

Suggested exploration: 80 Hz (large diagonal line near the first mode), 440 Hz (narrow ellipse), 700 Hz and 1000 Hz (lines and ellipses with different tilt), and 5–20 kHz (the pattern shrinks to a jittering dot). Then switch the waveform to Square or Sawtooth at about 30–150 Hz to see harmonics excite the mirror modes. Use **Kick** to watch a transient ring down into the steady pattern.

---

## 2. Coordinate system and geometry (Panel E, 3D tab)

* World axes: **x** to the right, **y** along the speaker axis (out of the cone, toward the laser), **z** up. Units are SI inside the code.
* Speaker centre: (0, 0, 0.165) m. Woofer cone radius is 75 mm. The tweeter has a 22 mm radius and sits at the centre.
* The **mirror** (radius 6 mm) sits at (−35, 0, 165) mm, **to the left of the tweeter**, facing +y.
* The **spring** (L = 25.4 mm) lies along +x on the cone surface. Its root is held by **blue tape** on the cone at x ≈ −63 mm. Its tip is **directly beneath** the mirror centre (offset e_v below it and e_h to its left). The mirror rests on the tip.
* The **laser** sits at (120, 260, 45) mm: the front-right corner of the wooden base, raised 4.5 cm and aimed at the mirror. The **speaker panel** stands perpendicular to the base.
* The **screen** is a plane perpendicular to the nominal reflected ray, at distance D from the mirror. Screen axes: u = (d_out0 × ẑ)/|·| is horizontal and v = u × d_out0 is vertical.

The 3D schematic uses the **same `G` geometry struct and the same `laserSpot` function** as the projection. The rays drawn there are therefore the computed rays. A "Tilt drawn ×10/×30" option exaggerates the drawing only; the projection always uses the true angles.

---

## 3. Governing equations

### 3.1 Tone generator (Panel A)
The band-limited Fourier sine series is

u(t) = Σₙ bₙ sin(2π n f t), for n f ≤ 22,050 Hz (a 44.1 kHz audio band). At most 300 harmonics are used.

| waveform | bₙ |
|---|---|
| sine | b₁ = 1 |
| square | 4/(πn), n odd |
| triangle | (8/π²n²)(−1)^((n−1)/2), n odd |
| sawtooth | (2/π)(−1)^(n+1)/n |

Notes use 12-tone equal temperament, f(n) = 440·2^((n−69)/12), with MIDI n = 21 (A0) … 119 (B8).

### 3.2 Amplifier (Panel B)
The signal path is volume × low-shelf × high-shelf (first-order shelving filters, s = iω):

* Bass: H_L(s) = (s + ω_L√G_L)/(s + ω_L/√G_L), with f_L = 100 Hz. The gain is G_L at DC and 1 at high frequency.
* Treble: H_H(s) = G_H (s + ω_H/√G_H)/(s + ω_H√G_H), with f_H = 10 kHz. The gain is 1 at DC and G_H at high frequency.
* Volume: linear amplitude gain V/100 (assumed taper).

For a pure sine, bass and treble change **only the amplitude and phase of that one frequency**. For square, triangle and sawtooth waves they also change the relative harmonic amplitudes. **They never create a second frequency.**

### 3.3 Speakers (Panel C)
Each cone is a driven damped oscillator, mᵢ s̈ᵢ + cᵢ ṡᵢ + kᵢ sᵢ = Fᵢ(t). For a harmonic drive this gives the displacement FRF

sᵢ = Aᵢ · Hᵢ(f) · U, with Hᵢ = 1 / (1 − r² + 2iζ_s r) and r = f/f_s,i.

Here Aᵢ is the static displacement at full drive. Both speakers receive **the same U(f)**. They differ in A, f_s, and an optional extra tweeter phase φ (polarity/wiring).

The woofer is the mechanical support. Its acceleration ÿ_b = −ω²s_out is the base excitation of the spring. The tweeter is **not** assumed to be rigidly attached to the spring. It acts on the spring through two modelled paths:

1. **Structural coupling κ.** The tweeter's reaction rocks the cone locally at the tape, by the angle ψ_b = κ s_in / r_tape (r_tape = 50 mm).
2. **Acoustic near field.** The tweeter acts as a baffled simple source, p = ρ₀ S_tw s̈_in / (2π d), with d = 35 mm. This pressure pushes on the mirror face, whose centre of pressure is shifted toward the tweeter.

### 3.4 Tape, spring and mirror (Panel D): 3-DOF base-excited model
The generalised coordinates are **q = [z, θ_v, θ_h]**:

* z is the spring-tip translation along y, relative to the cone (the "relative displacement" of the 1-DOF model).
* θ_v (displayed as **θ_y, vertical tilt**) is the rotation about the horizontal x axis, which is spring torsion. It moves the spot vertically.
* θ_h (displayed as **θ_x, horizontal tilt**) is the rotation about the vertical z axis, which is spring bending. It moves the spot horizontally.

The equation of motion is

**M q̈ + C q̇ + K q = −M r_y ÿ_b − M r_ψ ψ̈_b − r_p F_ac**

* **Mass matrix (exact small-angle rigid body).** The mirror centre of mass sits at offset (e_h, 0, e_v) from the spring tip. Its y-velocity is ż + e_h θ̇_h − e_v θ̇_v, so M = m b bᵀ + diag(0, I_c, I_c), with b = [1, −e_v, e_h]ᵀ and I_c = mR²/4. This gives J_v = I_c + m e_v² and J_h = I_c + m e_h² as rotational inertias, plus the translation–rotation coupling terms.
* **Spring stiffness.** The coil spring is treated as an Euler–Bernoulli cantilever. Its tip stiffness in bending is (EI/L³)[12, −6L; −6L, 4L²] acting on (z, θ_h), with **k = 3EI/L³** (the slider value). The torsion stiffness is k_θv = β·EI/L, with β = 0.3 (estimate). The cross-axis term is K₂₃ = κ_c √(k_θv · EI/L). It represents the slant of the cone surface and helix bend–twist coupling, and K stays positive-definite for |κ_c| < 1.
* **Tape / contact compliance.** The tape and the mirror-on-spring contact form a series compliance, K = (K_s⁻¹ + C_t)⁻¹, with C_t = diag(c, c/w², c/w²) and w = 10 mm. This is the simplified compliant-contact model: the joint is a linear series spring, and lift-off is handled by the contact limit below.
* **Damping.** Modal (proportional) damping with ratio ζ: C = MΦ diag(2ζωᵢ) ΦᵀM.
* **Influence vectors.** r_y = [1, 0, 0]ᵀ is cone translation. r_ψ = [L, 0, 1]ᵀ is cone rocking about z at the spring root. r_p = [1, −e_v, e_h + 0.3R]ᵀ is the acoustic force point. Each generalised force has the units of its coordinate: N for z and N·m for θ.
* **Reference 1-DOF model.** m z̈ + c ż + k_eff z = −m ÿ_b, with f_n = (1/2π)√(k_eff/m) and k_eff = (1/k + c_t)⁻¹. Statically condensing the full bending matrix (θ_h free) returns exactly k, so the 3-DOF model reduces to this one.
* **Contact limit (nonlinearity).** Before the optics, θ → θ_lim·tanh(θ/θ_lim), with θ_lim = 0.15 rad. A warning appears when the tilt approaches the limit, because the real mirror would lift or rattle there.

**Why lines and ellipses appear.** In steady state every coordinate oscillates at the drive frequency (and its harmonics). The figure is a **line** when θ_x and θ_y are in phase or anti-phase. It becomes an **ellipse** when damping near a mode, or the second (tweeter) excitation path with its different phase, shifts one axis relative to the other. The ellipse tilt follows the amplitude ratio. Rotating or evolving figures only come from real mechanisms in the model: decaying transients at the natural frequencies (Reset/Kick or any parameter change), harmonics of non-sine waveforms, the contact nonlinearity, and the random base disturbance.

### 3.5 Solution method (real time)
* **Steady state.** For each harmonic, Q_n = Φ[(ΦᵀF_n)/(ω_i² − ω² + 2iζω_iω)], which is exact. The result over one drive period is stored in a lookup table of 256–4096 points, so each frame is a cheap periodic interpolation.
* **Transient.** This is the exact modal free response, η_i(τ) = e^{−ζω_iτ}[η₀cos ω_dτ + (η̇₀ + ζω_iη₀)/ω_d sin ω_dτ]. Whenever a parameter changes at time t_c, the old state at t_c becomes the new initial condition, so the motion is continuous. The difference from the new steady state then rings down at the natural frequencies.
* **Disturbance.** Narrow-band random base acceleration (12 sinusoids clustered within ±ζω_i of each mode, fixed random phases) passed through the same linear response.
* **No ODE in the animation loop.** `ode45` is used only in verification test 11, to confirm the analytical solution.
* **Sampling and anti-aliasing.** Samples resolve 64 points per drive period, 12 per highest retained harmonic, and 24 per mirror mode. If the exposure window would need more than ~3000 samples (high frequencies), only the most recent part of the window is drawn. The pattern is periodic, so it looks the same, and it is never undersampled. The sample budget also adapts to the measured frame time.
* **Rendering.** A 30 ms `timer` advances simulated time by wall-time × speed. The sim speed is a slow-motion factor; at 0.01 the 440 Hz spot visibly travels the ellipse. Each frame only updates existing line and patch objects. Old parameter segments are kept for one simulated second so the trail stays continuous across slider moves.

### 3.6 Optics
* Mirror normal: n = R_z(θ_h) R_x(θ_v) n₀, the exact rotation.
* Mirror plane: through the displaced mirror centre (z + e_hθ_h − e_vθ_v along y).
* Hit point: the intersection of the laser ray with that plane.
* **d_out = d_in − 2(d_in·n)n**.
* Screen point: the intersection of the reflected ray with the screen plane, expressed as screen coordinates (u, v) in cm.
* Invalid geometry (beam hitting the mirror's back, or grazing the screen) gives NaN and is not drawn.

### 3.7 Idealized Lissajous model (Panel G, comparison only)
x = A sin(2π(a f + Δ)t + δ), y = B sin(2π b f t). The figure is labelled **IDEALIZED** on screen. The optional "bass(x)/treble(y)" mode scales A and B by the bass and treble gains and is labelled **HYPOTHETICAL two-frequency visualization**. It is *not* how an equaliser works and is never used by the physics model.

---

## 4. Controls

| Panel | Control | Meaning |
|---|---|---|
| A | Frequency (log slider), exact Hz field | 1–20,154 Hz. Both stay synchronised, and changes apply live. |
| A | ◀ / ▶, step | Decrease/increase by the step (0.001–100 Hz, default 1 Hz, like the website). |
| A | ×½ / ×2 | Down or up one octave. |
| A | Musical note | A0…B8 in 12-TET. Choosing a note sets f, and typing f shows the nearest note with cents. |
| A | Waveform | Sine, square, triangle or sawtooth (Fourier series above). |
| B | Master volume | 0–100 %, linear amplitude. |
| B | Bass / Treble | ±12 dB shelving EQ at 100 Hz / 10 kHz. The readout shows the EQ gain at f. |
| C | Outer / inner amplitude | Static cone displacement at full drive. |
| C | Speaker damping ζ_s | Damping of both cone oscillators. |
| C | Relative phase φ | Extra tweeter phase (polarity or wiring). The frequency is still the same. |
| C | Woofer / tweeter f_s | Cone resonance frequencies. |
| C | Speaker coupling κ | Tweeter-induced rocking of the tape mount. |
| D | Spring stiffness k | Tip bending stiffness 3EI/L³ (N/m). |
| D | Mirror mass m | Mirror + adhesive mass (inertia scales with it). |
| D | Damping ratio ζ | Modal damping of the spring–mirror modes. |
| D | Tape/contact compliance | Series compliance of tape and contact. Larger values soften every mode. |
| D | Horizontal / vertical tilt sensitivity | CoM lever arms e_h, e_v (mm) that turn cone acceleration into torque about z / x. |
| D | Cross-axis coupling κ_c | Stiffness coupling between the two tilts. |
| D | Initial perturbation θ₀ | Angular offset applied at Reset and on Kick. |
| D | Readouts and lamp | 1-DOF f_n, coupled mode frequencies, and resonance indicator (red: within 2ζ of a mode; yellow: within 10ζ). |
| E | Mirror–screen distance D | Moves the screen plane. Pattern size scales with D. |
| E | Brightness, spot size | Laser rendering. |
| E | Projection scale | Visible half-width of the screen in cm (used when auto-fit is off). |
| E | Exposure / persistence | Length of the visible trail (simulated seconds), like a camera exposure. |
| E | Disturbance | RMS of the random base vibration (m/s²). Set it to 0 for a perfectly stable figure. |
| E | Auto-fit | Zooms the view to the pattern. True sizes are always printed on screen. |
| F | Play / Pause / Reset / Kick / Verify | Animation control. Reset restarts the sound with the mirror at rest; Kick applies an angular tap. |
| F | Simulation speed | Slow-motion factor from 10⁻⁴ to 1. |
| F | Display model | Physics model or idealized Lissajous. |
| F | Update 3D schematic | Turns schematic refresh on or off (for performance). |
| G | a, b, δ, A, B, detune, hypothetical | Idealized Lissajous parameters. |

Right-hand tabs:
* **Time signals**: speaker displacements, spring-tip displacement z(t), and mirror tilts θ_x(t), θ_y(t) over the last 3 periods.
* **Resonance**: steady-state |θ_x|, |θ_y|, |z| and spot size vs drive frequency (log–log) at the current amplifier settings, with the modal frequencies and the current f marked.
* **3D schematic**: wooden base and panel, woofer cone, tweeter, blue tape, 1" coil spring, tilting mirror with its normal, laser pointer, incoming and reflected rays, and the projection screen with the spot.

The default constants give coupled modes at about **86, 433 and 952 Hz**, with a 1-DOF reference f_n ≈ 130 Hz.

---

## 5. Numerical verification (button **Verify**, or `LaserSpeakerSimulator('verify')`)

| # | Check | Pass criterion |
|---|---|---|
| 1 | k×2 raises f_n by exactly √2 (1-DOF and all coupled modes, no tape); k×1.5 with tape raises all modes | error < 1e-7 |
| 2 | m×2 lowers all natural frequencies by exactly 1/√2 | error < 1e-7 |
| 3 | 3-DOF resonant peak falls as ζ = 0.01 → 0.03 → 0.1; 1-DOF base-excited peak equals 1/(2ζ√(1−ζ²)) | monotone, rel. error < 1e-6 |
| 4 | Phase lag: −90° exactly at r = 1, monotone in r; 3-DOF z phase changes by more than 90° across the first mode | error < 1e-9 |
| 5 | 2000 random rays: \|d_out\| = 1, equal angles, coplanarity | < 1e-12 |
| 6 | Mirror rotation δ about the axis normal to the plane of incidence turns the beam by 2δ | \|ratio − 2\| < 1e-6 |
| 7 | Spot displacement at D = 1, 2, 4 m grows (≈ proportional) | ratio error < 1 % |
| 8 | Equal-frequency Lissajous: δ=0 line, δ=90° circle, δ=45° axis ratio tan 22.5°, A=2B ellipse | < 1e-6 |
| 9 | Rational ratios 3:2, 5:4, 1:3, 2:1 close after T = 1/gcd; √2:1 does not | < 1e-9 |
| 10 | 4 waveforms × 33 frequencies from 1 to 20,154 Hz: finite spots, bounded samples and table sizes, slider round trip | 0 invalid |
| 11 | Analytical steady state + transient vs `ode45` (60 Hz switched on from rest); lookup-table error | rel. error < 1e-4 |
| 12 | 12-TET note table (A0, C4, A4, B8) and nearest-note detection | < 1e-9 |

These tests show **internal numerical consistency** only. They do not show that the model matches the physical experiment.

> I could not run MATLAB in the environment where this was written, so the code has not been executed. The physics model and default constants were checked with an independent numerical prototype (mode frequencies, amplitudes and phases quoted above). The MATLAB source was reviewed by hand and checked for balanced blocks. Run `LaserSpeakerSimulator('verify')` first on your machine. If any API differs in your MATLAB release, the error message will name the line.

---

## 6. Assumptions and limitations

* **Uncalibrated constants.** These were not measured from the photographs: spring stiffness, torsion ratio, mirror mass, lever arms, tape compliance, speaker f_s / damping / excursion, coupling κ, acoustic distance and EQ corner frequencies. Treat every number as a placeholder.
* **Speakers.** Each is a single linear piston mode with no voice-coil inductance, crossover, cone break-up, or enclosure/air-load model. The cone motion at the tape is assumed equal to the piston motion plus one rocking term.
* **Spring.** A coil spring is idealised as an Euler–Bernoulli cantilever plus a torsion spring. Its axial mode, self-mass, and higher (surge) modes are neglected. The mirror is a rigid disc.
* **Contact.** It is a linear series compliance plus a smooth tanh angular limit. There is no impact or stick–slip model, and the tanh is applied at the output (it does not feed back into the dynamics).
* **Linear superposition** of steady state, transient and disturbance. Real large-amplitude behaviour (rattling, spring buckling, tape creep) is not modelled.
* **Acoustic path** uses a simple-source near-field formula. There is no diffraction, cavity or room acoustics.
* **Optics.** The beam is an ideal ray (no divergence or speckle), the mirror is flat and front-surface, and the screen is flat and perpendicular to the nominal ray. A camera's response (exposure, saturation, rolling shutter) is approximated only by the persistence window and the glow rendering.
* **Disturbance** is a stationary narrow-band random base vibration with fixed seeds, not a measured spectrum.

---

## 7. Comparing the simulation with experimental photographs

1. **Calibrate the geometry first.** Measure laser, mirror and screen positions and D. Enter them in `defaultParams` (`P.geo.*`, `P.D`) and the screen-plane orientation follows. Put a ruler on the screen so photo pixels convert to cm.
2. **Measure the pattern in each photo.** Fit an ellipse to the trace (for example, threshold the red channel and fit by least squares or with `regionprops` / principal axes). Record the major length, minor width, major-axis angle and axis ratio. The simulator prints the same quantities ("pattern L × W cm, major axis …°", plus the axis ratio in the status line), computed the same way from the trail.
3. **Match the camera exposure.** Set "Exposure" to the camera shutter time. At audio frequencies this makes little difference, but it matters below about 20 Hz.
4. **Sweep the frequency experimentally.** Record pattern size vs frequency, ideally with a fixed volume and a sine. Compare the curve with the **Resonance** tab's spot-size line. Peak frequencies give the natural frequencies, and fitting them sets k, β_t, m and C_tape. Peak widths (half-power bandwidth Δf ≈ 2ζf_n) give ζ.
5. **Measure the phase.** The minor/major ratio and tilt of an ellipse give the phase between θ_x and θ_y: for x = A sin(ωt + δ), y = B sin ωt, the ellipse area is πAB|sin δ|. Use that phase to set κ, φ and e_h / e_v.
6. **Measure the parts directly where you can.** Weigh the mirror. Find the spring stiffness by hanging small masses. Measure speaker excursion with a phone slow-motion video or an accelerometer app. Measure the free ring-down after a tap (the Kick button) to get the modal frequencies and decay rate.
7. **Check the waveforms.** At low frequency, square vs sine should change the trace shape the way the harmonic model predicts. If the measured trace does not change, the speaker or amplifier is probably band-limiting the harmonics.
8. **Report residuals,** not just visual similarity. Give the per-frequency error in size, aspect ratio and angle, and keep the parameters labelled as fitted estimates.
