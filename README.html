<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
@page { size: Letter; margin: 0.42in 0.45in 0.42in 0.45in; }
body { font-family: "DejaVu Sans", Arial, sans-serif; font-size: 7.55pt; line-height: 1.27; color: #111; margin: 0; }
h1 { font-size: 13pt; margin: 0 0 1px 0; color: #8b0000; }
.sub { font-size: 7.5pt; color: #444; margin-bottom: 4px; }
h2 { font-size: 9pt; margin: 6px 0 2px 0; color: #fff; background: #8b0000; padding: 1.5px 5px; }
table { border-collapse: collapse; width: 100%; margin: 0 0 2px 0; }
th, td { border: 0.6px solid #999; padding: 1.6px 4px; vertical-align: top; }
th { background: #eee; text-align: left; font-size: 7.4pt; }
td.ln, th.ln { width: 13%; font-family: "DejaVu Sans Mono", monospace; font-size: 7.1pt; }
td.eq { font-family: "DejaVu Serif", "Times New Roman", serif; font-size: 8.1pt; }
.m { font-family: "DejaVu Sans Mono", monospace; font-size: 7pt; }
.flow { border: 0.8px solid #8b0000; padding: 3px 6px; margin: 3px 0 2px 0; background: #fff6f6; font-size: 7.6pt; }
.note { font-size: 7pt; color: #333; }
.pb { page-break-before: always; }
i { font-family: "DejaVu Serif", serif; }
.two { display: flex; gap: 8px; }
.two > div { flex: 1; }
</style></head><body>

<h1>Sound-Driven Laser Pattern Simulator: Model Summary</h1>
<div class="sub">Line numbers refer to <span class="m">LaserSpeakerSimulator.m</span> exactly as delivered (2016 lines). All constants are illustrative estimates, not calibrated values. SI units are used inside the code.
Lines 1&ndash;1189: main function and GUI (nested functions). Lines 1196&ndash;2016: physics, optics and verification (local functions). If your numbering is offset, use Ctrl+F on the function name given.</div>

<div class="flow"><b>Signal chain.</b> Tone <i>u</i>(<i>t</i>) <b>></b> amplifier (volume, bass, treble) <b>></b> woofer + tweeter cones <i>s</i><sub>out</sub>, <i>s</i><sub>in</sub> (same input) <b>></b> forcing on the taped spring <b>></b> mirror state <b>q</b> = [<i>z</i>, &theta;<sub>v</sub>, &theta;<sub>h</sub>] (nonlinear ODE, RK4) <b>&gt;</b> mirror normal <b>n</b>(<i>t</i>) <b>></b> reflection <b>></b> screen point (<i>X</i>, <i>Y</i>) <b>></b> 40&nbsp;ms persistence trail (one laser figure).
Default parameters: <span class="m">defaultParams</span>, lines <b>1196&ndash;1265</b>.</div>

<h2>1. Tone generator (website-style controls)</h2>
<table>
<tr><th style="width:33%">Equation</th><th>Variables</th><th class="ln">Code lines</th></tr>
<tr><td class="eq"><i>u</i>(<i>t</i>) = &Sigma;<sub><i>n</i></sub> <i>b<sub>n</sub></i> sin(2&pi;<i>n</i>&phi;), &nbsp; d&phi;/d<i>t</i> = <i>f</i><br>
sine: <i>b</i><sub>1</sub>=1; square: 4/(&pi;<i>n</i>), <i>n</i> odd<br>
triangle: 8(&minus;1)<sup>(<i>n</i>&minus;1)/2</sup>/(&pi;<sup>2</sup><i>n</i><sup>2</sup>), <i>n</i> odd<br>
sawtooth: 2(&minus;1)<sup><i>n</i>+1</sup>/(&pi;<i>n</i>); &nbsp; <i>n f</i> &lt;= 22&nbsp;050&nbsp;Hz</td>
<td><i>f</i> drive frequency, 1&ndash;20&nbsp;154&nbsp;Hz (<span class="m">P.f</span>); &phi; tone phase in cycles (<span class="m">S.phase</span>), advanced by <i>f</i>&middot;d<i>t</i> each step so frequency changes never make the waveform jump; <i>b<sub>n</sub></i> Fourier coefficients (<span class="m">b</span>); at most 512 harmonics (<span class="m">P.nHarmMax</span>).</td>
<td class="ln">waveformHarmonics 1349&ndash;1369 (b<sub>n</sub>: 1358, 1361, 1363)<br>phase advance 1507<br>UI panel A 230&ndash;269</td></tr>
<tr><td class="eq"><i>f</i>(<i>n</i>) = 440&middot;2<sup>(<i>n</i>&minus;69)/12</sup>; &nbsp; cents = 100(<i>n</i><sub>exact</sub> &minus; round)</td>
<td>MIDI note <i>n</i> = 21 (A0) &hellip; 119 (B8); nearest note and cent offset shown for any typed frequency.</td>
<td class="ln">noteFrequency 1684&ndash;1686<br>nearestNote 1688&ndash;1693</td></tr>
<tr><td class="eq"><i>f</i> := <i>f</i>&middot;2<sup>&plusmn;<i>r</i>&Delta;<i>t</i>/12</sup></td>
<td>Automatic sweep: <i>r</i> in semitones per simulated second (<span class="m">P.sweepRate</span>); &Delta;<i>t</i> simulated time per frame.</td>
<td class="ln">onTick 713&ndash;743 (sweep: 722)</td></tr>
</table>

<h2>2. Stereo amplifier (equaliser)</h2>
<table>
<tr><td class="eq" style="width:33%"><i>U</i>(<i>f</i>) = (<i>V</i>/100)&middot;<i>H<sub>L</sub></i>&middot;<i>H<sub>H</sub></i>&middot;(&minus;<i>i b<sub>n</sub></i>)<br>
<i>H<sub>L</sub></i> = (<i>s</i>+&omega;<sub>L</sub>sqrt<i>G<sub>L</sub></i>)/(<i>s</i>+&omega;<sub>L</sub>/sqrt<i>G<sub>L</sub></i>)<br>
<i>H<sub>H</sub></i> = <i>G<sub>H</sub></i>(<i>s</i>+&omega;<sub>H</sub>/sqrt<i>G<sub>H</sub></i>)/(<i>s</i>+&omega;<sub>H</sub>sqrt<i>G<sub>H</sub></i>), &nbsp;<i>s</i>=2&pi;<i>if</i></td>
<td><i>V</i> master volume 0&ndash;100&nbsp;% (linear taper, <span class="m">P.volume</span>); <i>G<sub>L</sub></i>, <i>G<sub>H</sub></i> = 10<sup>dB/20</sup> bass and treble gains, &plusmn;12&nbsp;dB (<span class="m">P.bass_dB</span>, <span class="m">P.treble_dB</span>); &omega;<sub>L</sub> = 2&pi;&middot;100&nbsp;Hz, &omega;<sub>H</sub> = 2&pi;&middot;10&nbsp;kHz shelf corners. EQ only rescales frequencies already present; it never creates a second frequency.</td>
<td class="ln">eqResponse 1314&ndash;1321<br>drive phasor U 1396</td></tr>
</table>

<h2>3. Two concentric speakers (same input, different response) and distortion</h2>
<table>
<tr><td class="eq" style="width:33%"><i>m<sub>i</sub></i>s&#776;<sub><i>i</i></sub> + <i>c<sub>i</sub></i>s&#775;<sub><i>i</i></sub> + <i>k<sub>i</sub>s<sub>i</sub></i> = <i>F<sub>i</sub></i>(<i>t</i>) &nbsp;=><br>
<i>S<sub>i</sub></i> = <i>A<sub>i</sub> U</i> / (1 &minus; <i>r</i><sup>2</sup> + 2<i>i</i>&zeta;<sub>s</sub><i>r</i>), &nbsp;<i>r</i> = <i>f</i>/<i>f</i><sub>s,<i>i</i></sub><br>
<i>S</i><sub>in</sub> includes the factor <i>e</i><sup><i>i</i>&phi;<sub>rel</sub></sup></td>
<td><i>S</i><sub>out</sub>, <i>S</i><sub>in</sub> woofer and tweeter cone displacement phasors (m); <i>A</i><sub>out</sub> = 0.6&nbsp;mm, <i>A</i><sub>in</sub> = 0.15&nbsp;mm static displacement at full drive; <i>f</i><sub>s,out</sub> = 55&nbsp;Hz, <i>f</i><sub>s,in</sub> = 1800&nbsp;Hz cone resonances; &zeta;<sub>s</sub> = 0.4 speaker damping; &phi;<sub>rel</sub> extra tweeter phase (wiring/polarity).</td>
<td class="ln">speakerFRF 1309&ndash;1312<br>cone phasors 1397&ndash;1398<br>(resonance plot: 1343&ndash;1344)</td></tr>
<tr><td class="eq"><i>S</i>(2<i>f</i>) = <i>h</i><sub>2</sub><i>S</i><sub>1</sub><sup>2</sup>/|<i>S</i><sub>1</sub>|, &nbsp; <i>S</i>(3<i>f</i>) = <i>h</i><sub>3</sub><i>S</i><sub>1</sub><sup>3</sup>/|<i>S</i><sub>1</sub>|<sup>2</sup></td>
<td>Weak cone nonlinearity: <i>h</i><sub>2</sub> = 2&nbsp;%, <i>h</i><sub>3</sub> = 1&nbsp;% (<span class="m">P.hd2</span>, <span class="m">P.hd3</span>); <i>S</i><sub>1</sub> fundamental cone phasor. A lightly damped mirror mode (<i>Q</i> = 1/2&zeta; = 50) amplifies 2<i>f</i> or 3<i>f</i>, which gives 2:1 and 3:1 Lissajous-like figures.</td>
<td class="ln">addDistortion 1371&ndash;1390 (1387&ndash;1388)<br>called at 1399</td></tr>
</table>

<h2>4. Forcing on the taped spring (base excitation, tweeter coupling, acoustic pressure)</h2>
<table>
<tr><td class="eq" style="width:33%"><b>F</b> = &minus;<b>M r</b><sub>y</sub> <i>a</i><sub>out</sub> &minus; <b>M r</b><sub>&psi;</sub> &psi;&#776;<sub>b</sub> &minus; <b>r</b><sub>p</sub> <i>F</i><sub>ac</sub><br>
<i>a</i><sub>out</sub> = &minus;&omega;<sup>2</sup><i>S</i><sub>out</sub>, &nbsp; &psi;&#776;<sub>b</sub> = &kappa;(&minus;&omega;<sup>2</sup><i>S</i><sub>in</sub>)/<i>r</i><sub>tape</sub><br>
<i>p</i> = &rho;<sub>0</sub>&pi;<i>a</i><sub>tw</sub><sup>2</sup>(&minus;&omega;<sup>2</sup><i>S</i><sub>in</sub>)<i>e</i><sup>&minus;<i>i</i>&omega;<i>d</i>/<i>c</i></sup>/(2&pi;<i>d</i>)<br>
<i>F</i><sub>ac</sub> = &pi;<i>R</i><sup>2</sup><i>p</i></td>
<td><i>a</i><sub>out</sub> woofer acceleration at the tape (m/s&sup2;); &psi;<sub>b</sub> local cone rocking caused by the tweeter (rad), &kappa; = 0.6 coupling (<span class="m">P.kappa</span>), <i>r</i><sub>tape</sub> = 50&nbsp;mm; <i>p</i> tweeter near-field pressure (Pa), <i>a</i><sub>tw</sub> = 20&nbsp;mm tweeter radius, <i>d</i> = 35&nbsp;mm tweeter&ndash;mirror distance, <i>c</i> = 343&nbsp;m/s; <i>R</i> = 6&nbsp;mm mirror radius; influence vectors <b>r</b><sub>y</sub> = [1,0,0], <b>r</b><sub>&psi;</sub> = [<i>L</i>,0,1], <b>r</b><sub>p</sub> = [1, &minus;<i>e<sub>v</sub></i>, <i>e<sub>h</sub></i>+0.3<i>R</i>]. Force rows are N (for <i>z</i>) and N&middot;m (for &theta;).</td>
<td class="ln">forceFromCone 1323&ndash;1333<br>a: 1326&ndash;1327<br>&psi;&#776;: 1328<br>p, F<sub>ac</sub>: 1329&ndash;1331<br>F: 1332<br>vectors 1299&ndash;1301</td></tr>
</table>

<h2>5. Tape, 1-inch spring and mirror: structural matrices</h2>
<table>
<tr><td class="eq" style="width:33%"><b>M</b> = <i>m</i> <b>b b</b><sup>T</sup> + diag(0, <i>I<sub>c</sub></i>, <i>I<sub>c</sub></i>), &nbsp;<b>b</b> = [1, &minus;<i>e<sub>v</sub></i>, <i>e<sub>h</sub></i>], &nbsp;<i>I<sub>c</sub></i> = <i>mR</i><sup>2</sup>/4</td>
<td>Exact small-angle rigid-body mass matrix. <i>m</i> = 0.5&nbsp;g mirror + adhesive mass (<span class="m">P.m</span>); <i>e<sub>h</sub></i>, <i>e<sub>v</sub></i> = 2&nbsp;mm centre-of-mass offsets from the spring tip (horizontal and vertical tilt sensitivity, <span class="m">P.e_h</span>, <span class="m">P.e_v</span>).</td>
<td class="ln">buildStructure 1267&ndash;1307<br>M: 1273&ndash;1275</td></tr>
<tr><td class="eq"><b>K</b><sub>s</sub> = [4<i>k</i>, 0, &minus;2<i>kL</i>; 0, &beta;<sub>t</sub><i>kL</i><sup>2</sup>/3, <i>K</i><sub>23</sub>; &minus;2<i>kL</i>, <i>K</i><sub>23</sub>, 4<i>kL</i><sup>2</sup>/3]<br>
<i>K</i><sub>23</sub> = &kappa;<sub>c</sub>sqrt(<i>k</i><sub>&theta;v</sub><i>k</i><sub>&theta;h</sub>); &nbsp; <i>k</i> = 3<i>EI</i>/<i>L</i><sup>3</sup></td>
<td>Euler&ndash;Bernoulli cantilever tip matrix (bending: <i>z</i>, &theta;<sub>h</sub>) plus torsion (&theta;<sub>v</sub>). <i>k</i> = 400&nbsp;N/m tip stiffness (<span class="m">P.k</span>); <i>L</i> = 25.4&nbsp;mm spring length; &beta;<sub>t</sub> = 1.15 torsion/bending ratio (helical spring: 1+&nu;/2, <span class="m">P.beta_t</span>); &kappa;<sub>c</sub> = 0.2 cross-axis coupling (<span class="m">P.kappa_c</span>).</td>
<td class="ln">1277&ndash;1283</td></tr>
<tr><td class="eq"><b>K</b> = (<b>K</b><sub>s</sub><sup>&minus;1</sup> + <b>C</b><sub>t</sub>)<sup>&minus;1</sup>, &nbsp;<b>C</b><sub>t</sub> = diag(<i>c</i>, <i>c</i>/<i>w</i><sup>2</sup>, <i>c</i>/<i>w</i><sup>2</sup>)</td>
<td>Tape and mirror-on-spring contact as a series compliance: <i>c</i> = 3&nbsp;mm/N (<span class="m">P.C_tape</span>), <i>w</i> = 10&nbsp;mm lever arm. This puts the two tilt modes close together.</td>
<td class="ln">1284&ndash;1286</td></tr>
<tr><td class="eq"><b>K</b>&Phi; = <b>M</b>&Phi;&Omega;<sup>2</sup>, &nbsp; <b>C</b> = <b>M</b>&Phi; diag(2&zeta;&omega;<sub><i>j</i></sub>) &Phi;<sup>T</sup><b>M</b><br>
<i>f<sub>n</sub></i> = (1/2&pi;)sqrt(<i>k</i><sub>eff</sub>/<i>m</i>), &nbsp;<i>k</i><sub>eff</sub> = (1/<i>k</i> + <i>c</i>)<sup>&minus;1</sup></td>
<td>Modes &Phi; (mass-normalised) and natural frequencies &omega;<sub><i>j</i></sub> (defaults 86.6, 379, 418&nbsp;Hz); modal damping ratio &zeta; = 0.01 (<span class="m">P.zeta</span>); 1-DOF reference natural frequency <i>f<sub>n</sub></i> (about 96&nbsp;Hz).</td>
<td class="ln">eig 1288&ndash;1292<br>C: 1298<br>f<sub>n</sub>: 1302&ndash;1303</td></tr>
<tr><td class="eq"><i>k</i><sub>3</sub> = &eta;&middot;<i>k</i><sub>rock</sub>/&theta;<sub>ref</sub><sup>2</sup>, &nbsp;<i>k</i><sub>rock</sub> = (<i>K</i><sub>22</sub>+<i>K</i><sub>33</sub>)/2</td>
<td>Contact hardening coefficient: &eta; = 2&nbsp;% stiffness rise at &theta;<sub>ref</sub> = 10&nbsp;mrad (<span class="m">P.eta_nl</span>, <span class="m">P.theta_ref</span>).</td>
<td class="ln">1305&ndash;1306</td></tr>
</table>

<div class="pb"></div>

<h2>6. Equation of motion and real-time solver (the source of the rotating / "3D" figures)</h2>
<table>
<tr><th style="width:33%">Equation</th><th>Variables</th><th class="ln">Code lines</th></tr>
<tr><td class="eq"><b>M q&#776;</b> + <b>C q&#775;</b> + <b>K q</b> + <i>k</i><sub>3</sub>(&theta;<sub>v</sub><sup>2</sup>+&theta;<sub>h</sub><sup>2</sup>)[0, &theta;<sub>v</sub>, &theta;<sub>h</sub>]<sup>T</sup> = <b>F</b>(&phi;) &minus; <b>M r</b><sub>y</sub><i>a<sub>n</sub></i>(<i>t</i>)<br><br>
<b>q</b> = [<i>z</i>, &theta;<sub>v</sub>, &theta;<sub>h</sub>]<sup>T</sup></td>
<td><i>z</i> spring-tip displacement relative to the cone, along the speaker axis (m); &theta;<sub>v</sub> vertical tilt, rotation about the horizontal <i>x</i> axis (shown as &theta;<sub>y</sub>, rad); &theta;<sub>h</sub> horizontal tilt, rotation about the vertical <i>z</i> axis (shown as &theta;<sub>x</sub>, rad). State <span class="m">S.x</span> = [<b>q</b>; <b>q&#775;</b>]. The cubic term makes elliptical orbits precess at &Omega; = &beta;<i>ab</i>/(4&omega;), with &beta; = <i>k</i><sub>3</sub>/<i>J</i> and semi-axes <i>a</i>, <i>b</i>. It also bends the resonance peaks and causes jumps and up/down-sweep hysteresis.</td>
<td class="ln">integrateChunk 1455&ndash;1512<br>nonlinear term 1494&ndash;1495<br>accelerations 1496&ndash;1498</td></tr>
<tr><td class="eq">RK4: <b>x</b><sub><i>k</i>+1</sub> = <b>x</b><sub><i>k</i></sub> + (&Delta;<i>t</i>/6)(<b>k</b><sub>1</sub> + 2<b>k</b><sub>2</sub> + 2<b>k</b><sub>3</sub> + <b>k</b><sub>4</sub>)<br>
&Delta;<i>t</i> = min[1/(48<i>f</i>), 1/(16<i>n</i><sub>max</sub><i>f</i>), 1/(30<i>f</i><sub>fast</sub>)]<br>
<i>f</i><sub>fast</sub> = <i>f</i><sub>n,max</sub>sqrt(1 + 3&eta;(&theta;<sub>max</sub>/&theta;<sub>ref</sub>)<sup>2</sup>)</td>
<td>Fixed-step integration over only the newly elapsed simulated time each frame (wall time &times; <span class="m">P.speed</span>). <i>n</i><sub>max</sub> is the highest significant harmonic; <i>f</i><sub>fast</sub> is the amplitude-stiffened top mode. The step budget adapts to about 12&nbsp;ms per frame, and non-finite states are reset.</td>
<td class="ln">stage loop 1480&ndash;1504<br>update 1505&ndash;1508<br>stepSize 1447&ndash;1453<br>advance 745&ndash;782</td></tr>
<tr><td class="eq"><b>M</b><sup>&minus;1</sup><b>F</b>(&phi;) = one-period table (inverse FFT), linearly interpolated at phase &phi;</td>
<td>The forcing for all harmonics (fundamental, waveform harmonics, distortion) is precomputed once per parameter change, so the cost per step does not depend on the waveform. <span class="m">Mdl.A</span> = <b>M</b><sup>&minus;1</sup><b>K</b>, <span class="m">Mdl.B</span> = <b>M</b><sup>&minus;1</sup><b>C</b>.</td>
<td class="ln">buildModel 1392&ndash;1428 (table 1417&ndash;1422)<br>periodicTable 1430&ndash;1435<br>lookup 1485&ndash;1492</td></tr>
<tr><td class="eq"><i>a</i><sub><i>n</i>,<i>k</i>+1</sub> = &rho;<i>a</i><sub><i>n</i>,<i>k</i></sub> + &sigma;sqrt(1&minus;&rho;<sup>2</sup>)&xi;<sub><i>k</i></sub>, &nbsp; &rho; = <i>e</i><sup>&minus;2&pi;<i>f</i><sub>c</sub>&Delta;<i>t</i></sup></td>
<td>Random base vibration (Ornstein&ndash;Uhlenbeck process): &sigma; RMS acceleration, 0.01&nbsp;m/s&sup2; (<span class="m">P.noise</span>); <i>f</i><sub>c</sub> = 300&nbsp;Hz corner; &xi; standard normal samples.</td>
<td class="ln">761&ndash;762, 1478</td></tr>
</table>

<h2>7. Linear (small-signal) solutions: resonance plot, initial state, verification</h2>
<table>
<tr><td class="eq" style="width:33%"><b>Q</b>(&omega;) = &Phi; [&Phi;<sup>T</sup><b>F</b> / (&omega;<sub><i>j</i></sub><sup>2</sup> &minus; &omega;<sup>2</sup> + 2<i>i</i>&zeta;&omega;<sub><i>j</i></sub>&omega;)]<br>
&eta;<sub><i>j</i></sub>(&tau;) = <i>e</i><sup>&minus;&zeta;&omega;<sub><i>j</i></sub>&tau;</sup>[&eta;<sub>0</sub>cos&omega;<sub>d</sub>&tau; + (&eta;&#775;<sub>0</sub>+&zeta;&omega;<sub><i>j</i></sub>&eta;<sub>0</sub>)/&omega;<sub>d</sub> sin&omega;<sub>d</sub>&tau;]</td>
<td><b>Q</b> complex steady-state amplitudes of <b>q</b> (modal receptance); &eta;<sub><i>j</i></sub> exact free (transient) response of mode <i>j</i>, with &omega;<sub>d</sub> = &omega;<sub><i>j</i></sub>sqrt(1&minus;&zeta;<sup>2</sup>). Used for the Resonance tab, the Reset state and checking the integrator. The base-excited 1-DOF form is <i>m z&#776;</i> + <i>c z&#775;</i> + <i>kz</i> = &minus;<i>m y&#776;</i><sub>b</sub>.</td>
<td class="ln">modalSolve 1335&ndash;1339<br>harmonicResponse 1341&ndash;1347<br>ssStatePhase 1514&ndash;1520<br>transientModal 1534&ndash;1547<br>resonance plot 979&ndash;1009</td></tr>
</table>

<h2>8. Laser optics (vector reflection and screen intersection)</h2>
<table>
<tr><td class="eq" style="width:33%"><b>n</b> = <b>R</b><sub>z</sub>(&theta;<sub>h</sub>) <b>R</b><sub>x</sub>(&theta;<sub>v</sub>) <b>n</b><sub>0</sub><br>
<b>d</b><sub>out</sub> = <b>d</b><sub>in</sub> &minus; 2(<b>d</b><sub>in</sub>&middot;<b>n</b>)<b>n</b><br>
hit: <b>H</b> = <b>L</b><sub>0</sub> + <i>s</i><b>d</b><sub>in</sub>, &nbsp; <i>s</i> = (<b>P</b><sub>c</sub>&minus;<b>L</b><sub>0</sub>)&middot;<b>n</b> / (<b>d</b><sub>in</sub>&middot;<b>n</b>)<br>
screen: <b>P</b> = <b>H</b> + <i>t</i><b>d</b><sub>out</sub>, &nbsp; <i>t</i> = (<b>S</b><sub>0</sub>&minus;<b>H</b>)&middot;<b>w</b> / (<b>d</b><sub>out</sub>&middot;<b>w</b>)<br>
<i>X</i> = (<b>P</b>&minus;<b>S</b><sub>0</sub>)&middot;<b>u</b>, &nbsp; <i>Y</i> = (<b>P</b>&minus;<b>S</b><sub>0</sub>)&middot;<b>v</b></td>
<td><b>n</b><sub>0</sub> = [0,1,0] mirror normal at rest; <b>L</b><sub>0</sub> = (0.12, 0.26, 0.045)&nbsp;m laser position (front-right corner of the base, raised); <b>M</b><sub>0</sub> = (&minus;0.035, 0, 0.165)&nbsp;m mirror centre (left of the tweeter); <b>P</b><sub>c</sub> = <b>M</b><sub>0</sub> + (<i>z</i> + <i>e<sub>h</sub></i>&theta;<sub>h</sub> &minus; <i>e<sub>v</sub></i>&theta;<sub>v</sub>)<b>&#375;</b> moving mirror centre; <b>d</b><sub>in</sub> unit incoming beam direction; <b>w</b> screen normal (along the nominal reflected ray); <b>u</b>, <b>v</b> screen horizontal and vertical axes; <b>S</b><sub>0</sub> = <b>M</b><sub>0</sub> + <i>D</i><b>w</b>, with <i>D</i> = 2.5&nbsp;m mirror&ndash;screen distance (<span class="m">P.D</span>). Small-angle check: the beam turns 2&delta; for a mirror rotation &delta;, so spot displacement ~ 2<i>D</i>&delta;. Invalid geometry gives NaN and is not drawn; tilt is clamped at 0.3&nbsp;rad.</td>
<td class="ln">buildGeometry 1607&ndash;1619<br>rotateNormal 1626&ndash;1633<br>laserSpot 1635&ndash;1652<br>P<sub>c</sub>: 1638&ndash;1639<br>hit: 1640&ndash;1642<br>reflection: 1643 (also reflectDir 1621&ndash;1624)<br>screen: 1644&ndash;1651<br>qToScreen 1583&ndash;1590</td></tr>
</table>

<div class="two">
<div>
<h2>9. Display, persistence and pattern statistics</h2>
<table>
<tr><td>Ring buffer of every integrated sample (30&nbsp;000 samples); the trail is the last <i>T</i><sub>exp</sub> = 40&nbsp;ms (<span class="m">P.exposure</span>). Drawing is decimated without aliasing (>= 12 samples per period of the fastest content). Pattern length/width/angle come from the principal axes of the trail covariance: length = 2sqrt(2&lambda;<sub>1</sub>), width = 2sqrt(2&lambda;<sub>2</sub>).</td>
<td class="ln">bufAppend 802&ndash;813<br>bufRecent 815&ndash;828<br>renderFrame 838&ndash;953<br>selectSamples 1558&ndash;1574<br>patternStats 1666&ndash;1681</td></tr>
<tr><td>Idealized comparison (never used by the physics): <i>x</i> = <i>A</i> sin(2&pi;(<i>af</i>+&Delta;)<i>t</i>+&delta;), <i>y</i> = <i>B</i> sin(2&pi;<i>bft</i>).</td>
<td class="ln">lissajousTrail 1592&ndash;1605</td></tr>
<tr><td>3D schematic uses the same geometry and ray tracer.</td>
<td class="ln">1038&ndash;1151</td></tr>
</table>
</div>
<div>
<h2>10. Controls and verification</h2>
<table>
<tr><td>Kick: angular-velocity impulse &theta;<sub>0</sub>&omega;<sub>2</sub> (shows ringing and precession).</td><td class="ln">onKick 611&ndash;621</td></tr>
<tr><td>Parameter change: model rebuilt, state and tone phase kept.</td><td class="ln">setParam 533&ndash;558<br>processPending 560&ndash;576</td></tr>
<tr><td>Reset: start from the small-signal steady state.</td><td class="ln">resetSimulation 578&ndash;602</td></tr>
<tr><td>15 numerical checks (natural frequencies, damping, phase, reflection law, 2&times; deflection, screen distance, Lissajous shapes, full frequency range, RK4 vs exact and <span class="m">ode45</span>, notes, precession, distortion, phase continuity).</td><td class="ln">runVerification 1739&ndash;1979<br>precessionRate 1981&ndash;1999</td></tr>
</table>
</div>
</div>
<div class="note"><b>Limitations.</b> All constants are uncalibrated estimates. Distortion is a fixed fraction of the fundamental. Hardening is a smooth cubic term with no lift-off or impact. Not modelled: parametric (<i>f</i>/2) excitation, cone break-up, spring surge modes, air damping. Speakers are single piston modes; optics assume an ideal ray, a flat mirror and a flat screen. The verification checks internal consistency, not agreement with the experiment.</div>
</body></html>
