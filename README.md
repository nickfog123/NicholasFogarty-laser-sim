# Sound-Driven Laser Pattern Simulator: User Guide

This app simulates the speaker–spring–mirror laser experiment. You can install it on Windows with `MyAppInstaller.exe`, or run the MATLAB code `simulator.m`.

The chain it models is:

1. A tone drives the amplifier and both speakers.
2. The woofer cone shakes the taped 1-inch spring.
3. The mirror on the spring tilts.
4. The reflected laser draws one red figure on a simulated screen.

The figure is computed from the physics as the program runs. It changes with frequency, waveform, volume and the mechanical settings.

All mechanical values are estimates, not measurements of the real apparatus. The equations are in the separate model-summary PDF.

---

## 1. What you need

You can run the simulator in two ways:

| | Option A: Windows installer | Option B: MATLAB code |
|---|---|---|
| File | `MyAppInstaller.exe` | `simulator.m` |
| Needs MATLAB? | No | Yes, R2020b or newer (no toolboxes) |
| Operating system | Windows | Windows, macOS or Linux |
| Internet | Needed during installation | Not needed |
| Can edit the code or defaults? | No | Yes |

For either option, a screen of about 1600 × 950 pixels or larger works best. Smaller screens work if you maximize the window (see section 11).

## 2. Install

### Option A: Windows installer (no MATLAB needed)

1. Go to the release page: https://github.com/nickfog123/NicholasFogarty-laser-sim/releases/tag/v1.1.0
2. Under **Assets**, download **`MyAppInstaller.exe`**.
3. Double-click the downloaded file.
   - If Windows shows **"Windows protected your PC"**, click **More info → Run anyway**. Windows shows this for programs that aren't digitally signed.
   - If Windows asks for administrator permission, click **Yes**.
4. Follow the installer's steps. You can usually choose the install folder and whether to add a desktop shortcut.
5. The installer downloads and installs the free **MATLAB Runtime**, which the app needs instead of MATLAB. This is a large download (often a few gigabytes), so the step can take a while. Keep the computer online until it finishes.
6. When installation is complete, the app starts. Later, start it from the Start menu or the desktop shortcut.

To uninstall, go to **Windows Settings → Apps**, find the simulator and click **Uninstall**. The MATLAB Runtime is listed as its own entry. You can uninstall it too if no other program needs it.

### Option B: MATLAB code

1. Make a folder, for example `Documents\LaserSim`.
2. Put **`simulator.m`** in that folder. Either download it from the GitHub repository, or in MATLAB click **Home → New Script**, paste the code, and save it as `simulator.m`. The file name must be exactly `simulator.m`, all lowercase.
3. Set MATLAB's **Current Folder** to that folder (the folder bar above the Command Window), or run `addpath('full\path\to\LaserSim')`.

MATLAB runs the file by its file name, so the command is `simulator`. The editor may warn that the function name on line 1 doesn't match the file name. The code still runs. To remove the warning, change line 1 to `function out = simulator(varargin)`.

## 3. First run

### If you used the installer (Option A)

1. Start the app from the Start menu or the desktop shortcut. The first start can take 10–30 seconds while the MATLAB Runtime loads.
2. The window opens paused, showing the 440 Hz pattern. Press **▶ Play** to start the animation.
3. Optional: press **Verify** (Panel F) to run the 15 built-in numerical checks. A report window lists each check as PASS or FAIL.
4. Close the app with the window's close button.

### If you use MATLAB (Option B)

**Step 1: check the code.** Type this in the Command Window:

```matlab
R = simulator('verify');
```

This runs 15 numerical checks without opening the window. It usually takes under 20 seconds. Each check prints PASS or FAIL, and the last line should read `15 / 15 passed`. To see the results as a table, run `struct2table(R)`.

**Step 2: open the simulator.**

```matlab
simulator
```

The window opens paused, showing the 440 Hz pattern. Press **▶ Play** to start the animation.

To keep a handle to the window for saving pictures (section 10), start it with `app = simulator;` instead.

**Step 3: close it** with the window's close button. This also stops the animation timer. Close the window before starting a second copy.

The verification checks only test the code's internal consistency. They do not show that the model matches the real experiment.

---

## 4. The window

| Area | What it shows |
|---|---|
| Left column | Control panels **A–G**. Scroll down if they don't all fit. |
| Centre | The projection screen, in centimetres. The text at the top-left gives the frequency, note, waveform, model, and the size and angle of the figure. Warnings appear at the bottom-left. A status bar sits under the screen. |
| Right column | Three tabs: **Time signals**, **Resonance**, **3D schematic**. |

**How the laser is drawn.** The program calculates where the laser spot is at every instant. The trail is the spot's path over the last *Exposure* milliseconds (40 ms by default, roughly what your eye or a camera blends together). Older parts of the trail are darker, and the bright dot is where the spot is now. At normal speed (1×) the spot moves hundreds of times per second, so the dot jumps around the figure from frame to frame. That is expected. To follow the dot, use slow motion (section 5.5).

**Auto-fit** (Panel E, on by default) zooms the screen to fit the figure, so a 5 mm figure and a 50 cm figure both fill the view. The real size is always printed at the top-left. Turn Auto-fit off to see size changes directly.

---

## 5. Things to try

### 5.1 The default pattern
Press **▶ Play**. At 440 Hz you should see a small, narrow, tilted ellipse about a centimetre long. Drag the frequency slider and the figure changes while the animation runs.

### 5.2 The frequency sweep (turning, "3D-looking" figures)
1. Type **25** in **Exact frequency (Hz)** and press Enter.
2. Press **Sweep ▲**. The frequency rises by 1 semitone per second, or one octave every 12 s. Going from 25 to 700 Hz takes about a minute.
3. Press **Stop sweep** to hold any frequency, or **Sweep ▼** to go back down.

With the default settings, expect roughly the following:

| Frequency | What the laser does | Why |
|---|---|---|
| 25–80 Hz | A diagonal line rising to the right at about 25°, getting longer as the frequency rises. It is briefly longer near 28–31 Hz and 42–48 Hz. | The woofer shakes harder. Near those two frequencies, the speaker's 3rd and 2nd harmonics hit the 87 Hz mirror mode. They push the mirror along the same line, so the figure stays a line. |
| 85–100 Hz | The line grows very long (tens of cm), then suddenly collapses near 100 Hz. | This is the first mirror mode. The stiffening tape contact pushes the peak higher and makes it jump. |
| 100–175 Hz | The line shrinks to about 1 cm. | The frequency has passed the first mode. |
| 175–215 Hz | The line flattens a little and opens into small ellipses and 2:1 loops. | The 2nd harmonic reaches the 379 and 418 Hz tilt modes. |
| 220–340 Hz | A narrow ellipse whose long axis turns counter-clockwise, from about 40° to just past vertical. | The timing between the two tilt directions changes between the modes. |
| 350–400 Hz | A line leaning to the left that grows to about 15 cm near 395 Hz. | This is the 379 Hz tilt mode. |
| 410–450 Hz | The figure opens into an ellipse, then flips to a line falling to the right. | This is the 418 Hz tilt mode. |
| 450 Hz – 2 kHz | A short line, about 1 cm, slowly flattening toward horizontal. | The frequency is above all the mirror modes. |
| Above about 3 kHz | The figure shrinks toward a dot. | The mirror can't follow such fast shaking. |

These ranges come from a numerical test of the same equations. They shift if you change the Panel C or D settings, and the Resonance tab shows the natural frequencies for your current settings.

Sweeping down through 110 → 80 Hz, the big line may come back at a lower frequency than where it vanished on the way up. This happens because the stiffening contact makes the response depend on the sweep direction.

### 5.3 A 3:1 figure with a square wave
1. Set **Waveform** to **Square**.
2. Type **120** in Exact frequency, set **Sweep rate** to 0.1–0.2 st/s, and press **Sweep ▲**.
3. Between about 125 and 128 Hz, the line should open into a 3:1 figure (an S or N shape, or a three-lobed loop). The figure changes shape as the frequency moves, then closes back into a line. A weaker one may appear near 138–141 Hz.

This works because a square wave also contains 3 times its frequency. When that 3rd harmonic lands on a tilt mode, the mirror moves at f in one direction and at 3f in another direction. For a closer look, press **Stop sweep** inside that range and step through it using ◀ ▶ with the step set to 0.1 Hz.

### 5.4 Tap the mirror (Kick)
1. Set the frequency to 600 Hz, which is away from the resonances.
2. Set **Simulation speed** to 0.05×, **Damping ratio ζ** to 0.003, **Kick size θ₀** to about 25 mrad, and **Exposure** to 10 ms.
3. Press **Kick**.

The mirror rings in its two tilt modes. The ringing figure changes between a line and an ellipse and turns as it dies away, which takes a few seconds at 0.05×. Press Kick again to repeat. Afterwards, put the settings back or reopen the program.

### 5.5 Slow motion
Set **Simulation speed** to 0.01× and **Exposure** to about 5 ms. The dot now visibly draws the figure. At 0.001× it crawls.

### 5.6 Turn the nonlinear effects off
Set **2nd-harmonic distortion**, **3rd-harmonic distortion**, **Contact hardening η** and **Disturbance** all to 0. The model is then linear, and each frequency gives one fixed line or ellipse. This shows which behaviour comes from the nonlinear effects.

### 5.7 Idealized Lissajous comparison
Set **Display model** to **Idealized Lissajous (comparison)**. Panel G then draws x = A sin(2π(a·f + detune)·t + δ), y = B sin(2π·b·f·t). Things to try:
- a = b = 1 with δ = 0° (line), 90° (circle) or 45° (ellipse).
- a:b = 3:2 for a classic Lissajous figure.
- **x detune** = 0.5 Hz makes the figure turn continuously, one cycle every 2 s at 1× speed.

This is a pure mathematical picture, not the experiment. The **HYPOTHETICAL** box scales x by the bass gain and y by the treble gain, which is not how a real amplifier works. Set Display model back to **Physics** to return to the experiment.

---

## 6. Controls

Changes take effect while the animation runs. Moving any frequency control stops a running sweep.

### Panel A: Frequency and waveform
| Control | Range (default) | What it does |
|---|---|---|
| Frequency (log) | 1–20,154 Hz (440) | Logarithmic slider for the tone frequency. |
| Exact frequency (Hz) | 1–20,154 (440) | Type a value and press Enter. |
| ×½ ◀ ▶ ×2, step | step 0.001–100 Hz (1) | ×½ and ×2 halve or double the frequency. ◀ and ▶ move it down or up by one step. |
| Musical note (12-TET) | A0–B8 (A4) | Jumps to a note. It shows "(custom frequency)" when the frequency is between notes. The line underneath shows the nearest note and how far off it is, in cents. |
| Waveform | Sine, Square, Triangle, Sawtooth (Sine) | Shape of the tone. Square and sawtooth add strong harmonics. |

### Panel B: Stereo amplifier
| Control | Range (default) | What it does |
|---|---|---|
| Master volume | 0–100 % (50) | Overall drive. The figure size scales with it. |
| Bass (100 Hz shelf) | ±12 dB (0) | Boosts or cuts low frequencies. |
| Treble (10 kHz shelf) | ±12 dB (0) | Boosts or cuts high frequencies. |

Bass and treble only change the strength of frequencies already present. They never add a new frequency. The readout line shows the EQ gain at the current frequency.

### Panel C: Speakers (both get the same signal)
| Control | Range (default) | What it does |
|---|---|---|
| Outer (woofer) amp. | 0–2 mm (0.6) | How far the woofer cone moves at full volume below its resonance. This is the main source of shaking. |
| Inner (tweeter) amp. | 0–0.5 mm (0.15) | How far the tweeter cone moves. |
| Speaker damping ζ | 0.05–1.5 (0.40) | How sharp the speakers' own resonances are. |
| Relative phase φ | −180° to +180° (0) | Extra delay of the tweeter relative to the woofer (reversed wires = 180°). It changes ellipse shapes. |
| Woofer f_s / Tweeter f_s | 20–300 Hz (55) / 300–6000 Hz (1800) | The speakers' resonance frequencies. |
| Speaker coupling κ | 0–1.5 (0.6) | How much the tweeter rocks the cone where the tape is. |
| 2nd- / 3rd-harmonic distortion | 0–10 % (2 % / 1 %) | Extra cone motion at 2× and 3× the frequency. This causes the 2:1 figures. |

### Panel D: Tape, spring and mirror
| Control | Range (default) | What it does |
|---|---|---|
| Spring stiffness k | 20–5000 N/m (400) | A stiffer spring raises all natural frequencies. |
| Mirror mass m | 0.1–3 g (0.5) | A heavier mirror lowers all natural frequencies. |
| Damping ratio ζ | 0.002–0.5 (0.010) | Lower values give sharper resonances, longer ringing and more dramatic figures. |
| Tape/contact compl. | 0–5 mm/N (3) | How soft the tape and the mirror's seat are. Softer means lower tilt frequencies that sit closer together. |
| Horiz. / Vert. tilt sens. e_h, e_v | 0–10 mm (2 / 2) | How far the mirror's centre of mass sits from the spring tip. Larger values make the shaking tilt the mirror more in that direction. |
| Cross-axis coupling κ_c | −0.8 to +0.8 (+0.20) | How much one tilt direction drags the other along. It changes line angles and how open the ellipses are. |
| Torsion/bending β_t | 0.2–3 (1.15) | Twisting stiffness of the spring compared with its bending stiffness. Values near 1 put the two tilt modes close together. |
| Contact hardening η | 0–20 % at 10 mrad (2 %) | How much the contact stiffens as the mirror tilts. It causes the jumps, the turning after a kick, and the up/down-sweep differences. 0 makes the model linear. |
| Kick size θ₀ | 0–50 mrad (5) | Strength of the tap applied by the **Kick** button. |

The readout lines show the simple 1-DOF natural frequency, the three coupled natural frequencies, and the gap between the two tilt modes. The lamp is green away from a resonance, yellow within 10ζ of one (10 % at ζ = 0.01), and red within 2ζ.

### Panel E: Laser and screen
| Control | Range (default) | What it does |
|---|---|---|
| Mirror-screen dist. | 0.3–6 m (2.5) | Figure size scales with this distance. |
| Laser brightness / spot size | 10–100 % (90) / 1–8 px (3) | Appearance only. |
| Projection scale | ±0.05 to ±100 cm (±10) | Size of the visible screen. Only used when Auto-fit is off. |
| Exposure / persistence | 1–200 ms (40) | Length of the trail, like a camera's shutter time. |
| Disturbance (base) | 0–0.5 m/s² (0.010) | Random table or room vibration, which makes the figure shimmer. 0 gives a perfectly steady figure. |
| Auto-fit | on/off (on) | Zooms the view to the figure. |

### Panel F: Playback, sweep and simulation
| Control | What it does |
|---|---|
| ▶ Play / Pause | Starts or stops time. Changes made while paused take effect when you press Play. |
| Reset | Restarts from a steady state at t = 0 and stops any sweep. It keeps your current slider settings. Reopen the program to get the defaults back. |
| Kick | Taps the mirror (size set by Kick size θ₀ in Panel D). Starts playback if paused. |
| Verify | Runs the 15 checks and opens a report window, then resumes. |
| Sweep ▲ / Sweep ▼ / Stop sweep | Automatic smooth frequency sweep. It stops by itself at 1 Hz or 20,154 Hz. |
| Sweep rate | 0.1–12 semitones per second (default 1). It slows down along with Simulation speed. |
| Simulation speed | 0.001–1× (default 1×). 1× is real time, as your eye sees it; lower values give slow motion. |
| Display model | Physics model, or Idealized Lissajous (comparison only). |
| Update 3D schematic | Untick to make the program faster. |

### Panel G: Idealized Lissajous (only used when Display model is set to Lissajous)
| Control | Range (default) |
|---|---|
| x / y frequency ratio a, b | 1–9 (3 / 2) |
| Phase δ | 0–360° (90) |
| Amplitude A (x) / B (y) | 0.1–30 cm (8 / 8) |
| x detune | 0–5 Hz (0) |
| HYPOTHETICAL bass(x)/treble(y) | off |

---

## 7. Reading the screen

**Top-left text**
- Frequency, nearest note with its offset in cents, waveform, and sweep status.
- Which model is shown.
- `pattern L x W cm`: the full length and width of the trail.
- `major axis`: the angle of the figure's long axis, measured counter-clockwise from the screen's horizontal axis as drawn.
- `view ±…`: the half-width of the visible screen.

**Status bar (under the screen)**
- **Small-signal fundamental:** the predicted tilt amplitudes at the drive frequency, θx (horizontal tilt) and θy (vertical tilt) in mrad, and the phase difference between them. A phase near 0° or ±180° gives a line. ±90° with equal amplitudes gives a circle. Anything else gives an ellipse.
- **Trail:** the measured size, the axis ratio (0 = line, 1 = circle), the largest tilt, and the number of harmonics in the drive.

**Panel F readout:** simulated time, solver steps per frame, and the solver's time step.

**Warnings (bottom-left)**
| Warning | Meaning | What to do |
|---|---|---|
| CPU-limited: simulated time running at …× | The computer can't keep up with real time. The physics is still correct, just running slower. | Lower Simulation speed, untick Update 3D schematic, or ignore it. |
| High frequency: last … ms of the exposure drawn | The drawing limit was reached. The figure repeats, so it looks the same. | Nothing. |
| Tilt beyond 0.25 rad | The tilt is outside what the model handles. A real mirror would rattle. | Lower the volume or move away from the resonance. |
| Numerical safeguard: state restarted | Extreme settings made the calculation blow up, so it restarted. | Lower the volume or hardening, raise the damping, then press Reset. |
| Paused: the trail shows the motion before any change | Changes are applied when you press Play. | Press ▶ Play. |

**Right-hand tabs**
- **Time signals:** the last 20 ms (or 3 cycles, up to 0.5 s) of motion:
  - speaker cones in µm (woofer orange, tweeter blue)
  - spring-tip motion in µm (green)
  - mirror tilts in mrad (θx red, θy yellow)
- **Resonance:** how strongly the system responds to a pure sine at each frequency, for the current amplifier settings:
  - θx and θy in mrad, spring tip z in µm, spot half-size in cm (dashed). The units differ but share one log axis.
  - Dotted lines mark the natural frequencies, and the white line marks the current frequency.
  - Use this tab to find where the figures change. It is a linear calculation, so it does not show the jumps or the distortion effects.
- **3D schematic:** the wooden base, speakers, blue tape, spring, mirror (with its normal in green), the laser beams and the screen.
  - **View** switches between a close-up and the full beam path.
  - **Tilt drawn ×10 / ×30** exaggerates the mirror tilt in this drawing only.

---

## 8. Comparing with your experiment

1. Type the exact frequency you used on the tone generator.
2. Set **Mirror-screen dist.** to your measured distance.
3. Set **Exposure** to your camera's shutter time, which is listed in the photo's details (for example, 1/30 s = 33 ms).
4. Turn **Auto-fit** off and set **Projection scale** to half the width of the area shown in your photo, so the sizes compare directly.
5. Compare the length, width and angle in the top-left text with your photo. A photo taken from the other side of the screen is a mirror image, so compare the size of the angle and check its direction.
6. Adjust **Master volume** until the size matches at one frequency, then check other frequencies. If the frequencies where the figure changes don't match, adjust Panel D (spring, mass, tape, β_t) while watching the Resonance tab.
7. The simulated screen always faces the reflected beam. If your real screen is at an angle to the beam, the real figure is stretched in that direction.

---

## 9. Changing the start-up values (Option B only)

The installed app (Option A) always starts with the values it was built with, so use the sliders to change settings. To change the start-up values themselves, edit `simulator.m`.

All start-up values are in `defaultParams`, lines **1196–1265** of `simulator.m`. Edit them, save the file, and reopen the program.

Values in the code use SI units (kg, m, s, rad, and fractions instead of percentages), even where the slider shows g, mm, ms, mrad or %. Keep each value inside its slider's range.

| Setting | Line | Default in code |
|---|---|---|
| Frequency `P.f` | 1200 | `440` (Hz) |
| Waveform `P.wave` | 1201 | `'Sine'` |
| Volume `P.volume` | 1207 | `50` (%) |
| Distortion `P.hd2`, `P.hd3` | 1224–1225 | `0.02`, `0.01` |
| Spring stiffness `P.k` | 1228 | `400` (N/m) |
| Mirror mass `P.m` | 1230 | `0.5e-3` (kg) |
| Damping ratio `P.zeta` | 1232 | `0.01` |
| Tape compliance `P.C_tape` | 1233 | `3e-3` (m/N) |
| Contact hardening `P.eta_nl` | 1239 | `0.02` |
| Kick size `P.theta0` | 1241 | `5e-3` (rad) |
| Screen distance `P.D` | 1244 | `2.5` (m) |
| Exposure `P.exposure` | 1249 | `0.040` (s) |
| Disturbance `P.noise` | 1250 | `0.01` (m/s²) |
| Simulation speed `P.speed` | 1253 | `1` |
| Sweep rate `P.sweepRate` | 1254 | `1` (semitones/s) |

**Geometry** (lines 1261–1264, in metres). The origin is on the wooden base below the speaker centre, in the plane of the speaker rim. The axes are:
- x: horizontal, along the speaker panel (the laser is at +x, the mirror at −x)
- y: out of the speaker toward the laser
- z: up

| Line | Code | Meaning |
|---|---|---|
| 1261 | `P.geo.speakerCenter = [0; 0; 0.165];` | Speaker centre |
| 1262 | `P.geo.mirrorCenter = [-0.035; 0.0; 0.165];` | Mirror centre |
| 1263 | `P.geo.laserPos = [0.12; 0.26; 0.045];` | Laser aperture |
| 1264 | `P.geo.n0 = [0; 1; 0];` | Direction the mirror faces at rest |

The screen is always placed facing the reflected beam, at the mirror–screen distance.

To change the window size, edit `[30 30 1580 920]` on **line 124**. The numbers are left, bottom, width and height in pixels.

---

## 10. Saving pictures

**Installed app (Option A):** use the Windows Snipping Tool (**Windows + Shift + S**) for a screenshot, or **Windows + Alt + R** (Xbox Game Bar) to record a video.

**MATLAB (Option B):** start the program with a handle, then export the window:

```matlab
app = simulator;
% press Play, set up the figure you want, then:
exportapp(app.Figure, 'laser_440Hz.png')
```

You can also use your operating system's screenshot or screen-recording tool.

---

## 11. Troubleshooting

### Installer (Option A)
| Problem | Fix |
|---|---|
| "Windows protected your PC" | Click **More info → Run anyway**. Only do this for the file downloaded from the release page above. |
| The installer can't download the MATLAB Runtime | Check the internet connection and any firewall or proxy, then run `MyAppInstaller.exe` again. |
| Installation fails with a permissions error | Right-click `MyAppInstaller.exe` and choose **Run as administrator**, or choose an install folder you can write to. |
| The app takes a long time to open | The first start loads the MATLAB Runtime and can take 10–30 seconds or more. Later starts are usually faster. |
| The app doesn't open at all | Restart the computer and try again. If it still fails, uninstall and reinstall both the app and the MATLAB Runtime. |

### MATLAB (Option B)
| Problem | Fix |
|---|---|
| `Unrecognized function or variable 'simulator'` | Set the Current Folder to the folder containing `simulator.m`, and check the file name is exactly `simulator.m`. |
| The editor warns that the function name doesn't match the file name | The code still runs. To remove the warning, change line 1 to `function out = simulator(varargin)`. |
| An error mentions `WordWrap`, `Scrollable`, `uigridlayout` or `exportapp` | Your MATLAB is older than R2020b. |
| A timer keeps running after a crash | Run `delete(timerfind('Name','LaserSpeakerSimulatorTimer'))`, then close the window. |

### Both options
| Problem | Fix |
|---|---|
| The window doesn't fit on the screen | Maximize it. The centre screen shrinks first, and the left panels scroll. In MATLAB you can also edit line 124. |
| Nothing moves | Press **▶ Play**. Check that Master volume is above 0. |
| The figure is noisy or shimmering | Set Disturbance to 0. |
| It's only a line | That's expected at many frequencies. Try the sweep (5.2) or the square wave (5.3). |
| A popup says "Animation stopped: …" | Write down the message. Press Reset, then Play. If it repeats, report the message. |
| Verify shows FAIL | Write down the check number and its details. |
| Running slowly | Untick Update 3D schematic, lower Simulation speed, or use a sine wave. Square and sawtooth waves at high frequencies take the most computation. |

---

## 12. Limits

- All mechanical constants are uncalibrated estimates. The frequencies in section 5 move when you change them.
- Speaker distortion is a fixed percentage. Real distortion grows with volume.
- The contact stiffening is smooth. Real rattling or lift-off of the mirror is not modelled.
- The laser is treated as an ideal ray, the mirror and screen are flat, and the screen faces the beam.
