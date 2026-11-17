function out = LaserSpeakerSimulator(varargin)
%LASERSPEAKERSIMULATOR  Interactive sound-driven laser pattern simulator.
%
%   Simulates a two-speaker (coaxial woofer + central tweeter) apparatus in
%   which a small mirror rests on a ~1 inch (2.54 cm) coil spring that is
%   taped to the woofer cone. A red laser reflects off the mirror onto a
%   flat screen. Sound from a tone generator -> stereo amplifier (volume,
%   bass, treble) -> both speakers -> taped spring -> mirror tilt ->
%   vector reflection -> ONE continuously animated laser trajectory.
%
%   HOW TO RUN
%     1. Put this file in a folder on the MATLAB path (or cd into it).
%     2. Type   LaserSpeakerSimulator          (opens the GUI, press Play)
%     3. Headless numerical checks:  R = LaserSpeakerSimulator('verify');
%   Requires MATLAB R2020b or newer (uifigure / uigridlayout). No toolboxes.
%
%   MODEL SUMMARY (all constants are ESTIMATES unless measured)
%     Tone      : phase-continuous band-limited Fourier series of sine /
%                 square / triangle / sawtooth (frequency changes never
%                 make the waveform jump, as in a real tone generator).
%     Amplifier : first-order low-shelf (bass, 100 Hz) and high-shelf (treble,
%                 10 kHz) equaliser, linear master volume.
%     Speakers  : each cone is a driven damped oscillator
%                 m_i s_i'' + c_i s_i' + k_i s_i = F_i(t)
%                 plus small 2nd/3rd-harmonic distortion of the cone motion
%                 (real drivers and amplifiers are not perfectly linear).
%     Mount     : base translation of the woofer cone, local cone rocking
%                 caused by the tweeter, and the tweeter's near-field
%                 pressure on the mirror face (with acoustic delay).
%     Mirror    : 3-DOF rigid body on the cantilevered coil spring and a
%                 compliant tape/contact joint, q = [z; theta_v; theta_h]:
%                 M q'' + C q' + K q + f_nl(q) = F(t)
%                 f_nl = k3*(theta_v^2+theta_h^2)*[0; theta_v; theta_h]
%                 is an isotropic hardening of the compliant contact. It
%                 makes elliptical orbits precess and lets response
%                 ellipses open, tilt and jump as the frequency is swept.
%     Solution  : the full nonlinear model is integrated in real time with
%                 fixed-step RK4, only over the newly elapsed simulated
%                 time each frame (never re-integrating history). Linear
%                 frequency-domain solutions are used for the resonance
%                 plot, for initial conditions and for verification.
%     Optics    : exact vector law of reflection d_out = d_in - 2(d_in.n)n,
%                 ray/mirror-plane and ray/screen-plane intersections.
%
%   See README.md for the derivations, assumptions and limitations.

% ---------------------------------------------------------------------
% Headless verification mode
% ---------------------------------------------------------------------
if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1})) && strcmpi(varargin{1}, 'verify')
    R = runVerification(defaultParams());
    printVerification(R);
    if nargout > 0, out = R; end
    return
end

% ---------------------------------------------------------------------
% Shared application state (visible to all nested functions)
% ---------------------------------------------------------------------
P   = defaultParams();          % user parameters (SI units internally)
St  = buildStructure(P);        % spring-mirror structural matrices & modes
Mdl = buildModel(P, St);        % excitation tables for the current tone
G   = buildGeometry(P);         % optical geometry
S   = struct();                 % runtime state
S.t = 0;                        % simulated time (s)
S.phase = 0;                    % tone phase in cycles, continuous across frequency changes
S.x = zeros(6, 1);              % state [z; theta_v; theta_h; velocities]
S.an = 0;                       % current random base acceleration (m/s^2)
S.rng = RandStream('mt19937ar', 'Seed', 2026);
S.playing = false;
S.needRebuild = false;
S.geoDirty = false;
S.resDirty = true;
S.schemDirty = true;
S.freqWidgetsDirty = false;
S.frame = 0;
S.NsMax = 3000;                 % adaptive cap on drawn trail samples
S.stepBudget = 6000;            % adaptive cap on RK4 steps per frame
S.simDebt = 0;                  % requested simulated time not yet integrated
S.cpuLimited = false;
S.effSpeed = 1;
S.lastSteps = 0;
S.lastDt = 0;
S.ampEst = 0;                   % recent max tilt amplitude (rad), sets the step size
S.nanReset = false;
S.sweepDir = 0;                 % -1, 0, +1
S.lastTic = tic;
S.viewHalf = P.viewHalf_cm;
S.saturated = false;
S.Nb = 30000;                   % ring-buffer length (samples)
S.buf = nan(6, S.Nb);           % rows: t, z, theta_v, theta_h, s_out, s_in
S.bufHead = 0;
S.bufCount = 0;
S.timer = [];
S.verifyFig = [];
H   = struct();                 % graphics handles
CLR = themeColors();
GL  = greekChars();

buildUI();
resetSimulation();
updateReadouts();
updateResonancePlot();
drawSchematicStatic();
renderFrame(true);

S.timer = timer('Name', 'LaserSpeakerSimulatorTimer', ...
    'ExecutionMode', 'fixedSpacing', 'Period', 0.030, 'BusyMode', 'drop', ...
    'TimerFcn', @(~, ~) onTick(), 'ErrorFcn', @(~, ~) stopPlayback());

if nargout > 0
    out = struct('Figure', H.fig);
end

% =====================================================================
%                         NESTED FUNCTIONS (GUI)
% =====================================================================

    % -----------------------------------------------------------------
    % UI construction
    % -----------------------------------------------------------------
    function buildUI()
        H.fig = uifigure('Name', 'Sound-Driven Laser Pattern Simulator', ...
            'Position', [30 30 1580 920], 'Color', CLR.bg, ...
            'CloseRequestFcn', @(~, ~) onClose());

        mainG = uigridlayout(H.fig, [1 3]);
        mainG.ColumnWidth = {410, '1x', 490};
        mainG.RowHeight = {'1x'};
        mainG.BackgroundColor = CLR.bg;
        mainG.Padding = [6 6 6 6];
        mainG.ColumnSpacing = 6;

        % ---------------- left column: control panels ----------------
        left = uigridlayout(mainG, [7 1]);
        left.Layout.Row = 1; left.Layout.Column = 1;
        left.Scrollable = 'on';
        left.BackgroundColor = CLR.bg;
        left.RowSpacing = 6;
        left.Padding = [0 0 6 0];
        nRows = [8 4 10 13 7 8 7];
        left.RowHeight = num2cell(panelHeight(nRows));

        buildPanelA(makePanel(left, 1, 'A  Frequency & waveform (tone generator)', nRows(1)));
        buildPanelB(makePanel(left, 2, 'B  Stereo amplifier', nRows(2)));
        buildPanelC(makePanel(left, 3, 'C  Speaker dynamics (both speakers get the same signal)', nRows(3)));
        buildPanelD(makePanel(left, 4, 'D  Tape, 1-inch spring & mirror  [estimated constants]', nRows(4)));
        buildPanelE(makePanel(left, 5, 'E  Laser & screen geometry', nRows(5)));
        buildPanelF(makePanel(left, 6, 'F  Playback, sweep & simulation', nRows(6)));
        buildPanelG(makePanel(left, 7, 'G  Idealized Lissajous model (comparison only)', nRows(7)));

        % ---------------- centre: projection screen ------------------
        centre = uigridlayout(mainG, [2 1]);
        centre.Layout.Row = 1; centre.Layout.Column = 2;
        centre.RowHeight = {'1x', 46};
        centre.BackgroundColor = CLR.bg;
        centre.Padding = [0 0 0 0];
        buildProjectionAxes(centre);

        H.status = uilabel(centre, 'Text', '', 'FontColor', CLR.text, ...
            'FontName', 'Consolas', 'FontSize', 12, 'WordWrap', 'on', ...
            'BackgroundColor', CLR.panel);
        H.status.Layout.Row = 2; H.status.Layout.Column = 1;

        % ---------------- right column: scientific plots -------------
        tg = uitabgroup(mainG);
        tg.Layout.Row = 1; tg.Layout.Column = 3;
        buildSignalsTab(uitab(tg, 'Title', 'Time signals', 'BackgroundColor', CLR.bg));
        buildResonanceTab(uitab(tg, 'Title', 'Resonance', 'BackgroundColor', CLR.bg));
        buildSchematicTab(uitab(tg, 'Title', '3D schematic', 'BackgroundColor', CLR.bg));
    end

    function h = panelHeight(n)
        h = 26 .* n + 36;
    end

    function g = makePanel(parent, row, titleStr, nRows)
        p = uipanel(parent, 'Title', titleStr, 'BackgroundColor', CLR.panel, ...
            'ForegroundColor', CLR.accent, 'FontWeight', 'bold', 'FontSize', 12);
        p.Layout.Row = row; p.Layout.Column = 1;
        g = uigridlayout(p, [nRows 3]);
        g.ColumnWidth = {132, '1x', 78};
        g.RowHeight = repmat({22}, 1, nRows);
        g.BackgroundColor = CLR.panel;
        g.RowSpacing = 4;
        g.ColumnSpacing = 6;
        g.Padding = [6 4 6 4];
    end

    function lb = addLabel(g, row, col, txt, color)
        if nargin < 5, color = CLR.text; end
        lb = uilabel(g, 'Text', txt, 'FontColor', color, 'FontSize', 11);
        lb.Layout.Row = row; lb.Layout.Column = col;
    end

    function [sl, vl] = addSlider(g, row, labelStr, lims, val, fmtFcn, applyFcn)
        % Label | slider | value readout. applyFcn(v) receives the raw
        % slider value; fmtFcn(v) returns the readout text.
        addLabel(g, row, 1, labelStr);
        sl = uislider(g, 'Limits', lims, 'Value', min(max(val, lims(1)), lims(2)), ...
            'MajorTicks', [], 'MinorTicks', [], 'FontColor', CLR.text);
        sl.Layout.Row = row; sl.Layout.Column = 2;
        vl = uilabel(g, 'Text', fmtFcn(sl.Value), 'FontColor', CLR.value, ...
            'FontSize', 11, 'HorizontalAlignment', 'right');
        vl.Layout.Row = row; vl.Layout.Column = 3;
        sl.ValueChangingFcn = @(~, e) onSlider(e.Value, vl, fmtFcn, applyFcn);
        sl.ValueChangedFcn  = @(src, ~) onSlider(src.Value, vl, fmtFcn, applyFcn);
    end

    function onSlider(v, vl, fmtFcn, applyFcn)
        vl.Text = fmtFcn(v);
        applyFcn(v);
    end

    function b = addButton(g, row, col, txt, cb)
        b = uibutton(g, 'Text', txt, 'BackgroundColor', CLR.button, ...
            'FontColor', CLR.text, 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) cb());
        b.Layout.Row = row; b.Layout.Column = col;
    end

    function dd = addDropdown(g, row, labelStr, items, itemsData, val, cb)
        addLabel(g, row, 1, labelStr);
        dd = uidropdown(g, 'Items', items, 'ItemsData', itemsData, 'Value', val, ...
            'BackgroundColor', CLR.field, 'FontColor', CLR.text, ...
            'ValueChangedFcn', @(src, ~) cb(src.Value));
        dd.Layout.Row = row; dd.Layout.Column = [2 3];
    end

    % ---------------- Panel A: frequency & waveform --------------------
    function buildPanelA(g)
        [H.freqSlider, H.freqValLbl] = addSlider(g, 1, 'Frequency (log)', ...
            [0 log10(P.fMax)], log10(P.f), @(v) sprintf('%s Hz', formatFreq(10^v)), ...
            @(v) setFrequency(10^v, 'slider'));

        addLabel(g, 2, 1, 'Exact frequency (Hz)');
        H.freqEdit = uieditfield(g, 'numeric', 'Limits', [P.fMin P.fMax], ...
            'Value', P.f, 'ValueDisplayFormat', '%.3f', ...
            'BackgroundColor', CLR.field, 'FontColor', CLR.text, ...
            'ValueChangedFcn', @(src, ~) setFrequency(src.Value, 'edit'));
        H.freqEdit.Layout.Row = 2; H.freqEdit.Layout.Column = [2 3];

        bg = uigridlayout(g, [1 6]);
        bg.Layout.Row = 3; bg.Layout.Column = [1 3];
        bg.ColumnWidth = {'1x', '1x', '1x', '1x', 50, 70};
        bg.Padding = [0 0 0 0]; bg.ColumnSpacing = 4; bg.BackgroundColor = CLR.panel;
        addButton(bg, 1, 1, [GL.times char(189)], @() setFrequency(P.f / 2, 'button'));
        addButton(bg, 1, 2, char(9664), @() setFrequency(P.f - P.fStep, 'button'));
        addButton(bg, 1, 3, char(9654), @() setFrequency(P.f + P.fStep, 'button'));
        addButton(bg, 1, 4, [GL.times '2'], @() setFrequency(P.f * 2, 'button'));
        lb = uilabel(bg, 'Text', 'step', 'FontColor', CLR.text, 'HorizontalAlignment', 'right');
        lb.Layout.Row = 1; lb.Layout.Column = 5;
        H.stepDrop = uidropdown(bg, 'Items', {'0.001', '0.01', '0.1', '1', '10', '100'}, ...
            'ItemsData', [0.001 0.01 0.1 1 10 100], 'Value', 1, ...
            'BackgroundColor', CLR.field, 'FontColor', CLR.text, ...
            'ValueChangedFcn', @(src, ~) setStep(src.Value));
        H.stepDrop.Layout.Row = 1; H.stepDrop.Layout.Column = 6;

        [items, data] = noteItems();
        H.noteDrop = addDropdown(g, 4, 'Musical note (12-TET)', items, data, 69, ...
            @(m) onNoteSelected(m));
        H.noteInfo = addLabel(g, 5, [1 3], '', CLR.value);
        H.waveDrop = addDropdown(g, 6, 'Waveform', {'Sine', 'Square', 'Triangle', 'Sawtooth'}, ...
            {'Sine', 'Square', 'Triangle', 'Sawtooth'}, P.wave, @(w) setParam('wave', w, 'dyn'));
        H.freqBig = uilabel(g, 'Text', '', 'FontColor', CLR.accent, 'FontSize', 18, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        H.freqBig.Layout.Row = [7 8]; H.freqBig.Layout.Column = [1 3];
        refreshFrequencyWidgets('init');
    end

    % ---------------- Panel B: amplifier ------------------------------
    function buildPanelB(g)
        addSlider(g, 1, 'Master volume', [0 100], P.volume, @(v) sprintf('%.0f %%', v), ...
            @(v) setParam('volume', v, 'dyn'));
        addSlider(g, 2, 'Bass (100 Hz shelf)', [-12 12], P.bass_dB, @(v) sprintf('%+.1f dB', v), ...
            @(v) setParam('bass_dB', v, 'dyn'));
        addSlider(g, 3, 'Treble (10 kHz shelf)', [-12 12], P.treble_dB, @(v) sprintf('%+.1f dB', v), ...
            @(v) setParam('treble_dB', v, 'dyn'));
        H.eqInfo = addLabel(g, 4, [1 3], '', CLR.value);
    end

    % ---------------- Panel C: speakers -------------------------------
    function buildPanelC(g)
        addSlider(g, 1, 'Outer (woofer) amp.', [0 2], P.A_out * 1e3, @(v) sprintf('%.2f mm', v), ...
            @(v) setParam('A_out', v * 1e-3, 'dyn'));
        addSlider(g, 2, 'Inner (tweeter) amp.', [0 0.5], P.A_in * 1e3, @(v) sprintf('%.3f mm', v), ...
            @(v) setParam('A_in', v * 1e-3, 'dyn'));
        addSlider(g, 3, ['Speaker damping ' GL.zeta], [0.05 1.5], P.zeta_spk, @(v) sprintf('%.2f', v), ...
            @(v) setParam('zeta_spk', v, 'dyn'));
        addSlider(g, 4, ['Relative phase ' GL.phi], [-180 180], P.phi_rel_deg, @(v) sprintf('%+.0f%s', v, GL.deg), ...
            @(v) setParam('phi_rel_deg', v, 'dyn'));
        addSlider(g, 5, 'Woofer f_s', [20 300], P.fs_out, @(v) sprintf('%.0f Hz', v), ...
            @(v) setParam('fs_out', v, 'dyn'));
        addSlider(g, 6, 'Tweeter f_s', [300 6000], P.fs_in, @(v) sprintf('%.0f Hz', v), ...
            @(v) setParam('fs_in', v, 'dyn'));
        addSlider(g, 7, ['Speaker coupling ' GL.kappa], [0 1.5], P.kappa, @(v) sprintf('%.2f', v), ...
            @(v) setParam('kappa', v, 'dyn'));
        addSlider(g, 8, '2nd-harmonic distortion', [0 10], P.hd2 * 100, @(v) sprintf('%.1f %%', v), ...
            @(v) setParam('hd2', v / 100, 'dyn'));
        addSlider(g, 9, '3rd-harmonic distortion', [0 10], P.hd3 * 100, @(v) sprintf('%.1f %%', v), ...
            @(v) setParam('hd3', v / 100, 'dyn'));
        H.spkInfo = addLabel(g, 10, [1 3], '', CLR.value);
    end

    % ---------------- Panel D: spring & mirror ------------------------
    function buildPanelD(g)
        addSlider(g, 1, 'Spring stiffness k', [log10(20) log10(5000)], log10(P.k), ...
            @(v) sprintf('%.0f N/m', 10^v), @(v) setParam('k', 10^v, 'dyn'));
        addSlider(g, 2, 'Mirror mass m', [0.1 3], P.m * 1e3, @(v) sprintf('%.2f g', v), ...
            @(v) setParam('m', v * 1e-3, 'dyn'));
        addSlider(g, 3, ['Damping ratio ' GL.zeta], [log10(0.002) log10(0.5)], log10(P.zeta), ...
            @(v) sprintf('%.3f', 10^v), @(v) setParam('zeta', 10^v, 'dyn'));
        addSlider(g, 4, 'Tape/contact compl.', [0 5], P.C_tape * 1e3, @(v) sprintf('%.2f mm/N', v), ...
            @(v) setParam('C_tape', v * 1e-3, 'dyn'));
        addSlider(g, 5, 'Horiz. tilt sens. e_h', [0 10], P.e_h * 1e3, @(v) sprintf('%.1f mm', v), ...
            @(v) setParam('e_h', v * 1e-3, 'dyn'));
        addSlider(g, 6, 'Vert. tilt sens. e_v', [0 10], P.e_v * 1e3, @(v) sprintf('%.1f mm', v), ...
            @(v) setParam('e_v', v * 1e-3, 'dyn'));
        addSlider(g, 7, ['Cross-axis coupling ' GL.kappa '_c'], [-0.8 0.8], P.kappa_c, @(v) sprintf('%+.2f', v), ...
            @(v) setParam('kappa_c', v, 'dyn'));
        addSlider(g, 8, ['Torsion/bending ' GL.beta '_t'], [0.2 3], P.beta_t, @(v) sprintf('%.2f', v), ...
            @(v) setParam('beta_t', v, 'dyn'));
        addSlider(g, 9, ['Contact hardening ' GL.eta], [0 0.2], P.eta_nl, ...
            @(v) sprintf('%.0f%% @10mrad', 100 * v), @(v) setParam('eta_nl', v, 'dyn'));
        addSlider(g, 10, ['Kick size ' GL.theta '_0'], [0 50], P.theta0 * 1e3, @(v) sprintf('%.1f mrad', v), ...
            @(v) setParam('theta0', v * 1e-3, 'none'));
        H.fnInfo1 = addLabel(g, 11, [1 3], '', CLR.value);
        H.fnInfo2 = addLabel(g, 12, [1 3], '', CLR.value);
        lg = uigridlayout(g, [1 2]);
        lg.Layout.Row = 13; lg.Layout.Column = [1 3];
        lg.ColumnWidth = {22, '1x'}; lg.Padding = [0 0 0 0]; lg.BackgroundColor = CLR.panel;
        H.resLamp = uilamp(lg, 'Color', [0.2 0.8 0.3]);
        H.resLamp.Layout.Row = 1; H.resLamp.Layout.Column = 1;
        H.resText = uilabel(lg, 'Text', '', 'FontColor', CLR.text, 'FontSize', 11);
        H.resText.Layout.Row = 1; H.resText.Layout.Column = 2;
    end

    % ---------------- Panel E: optics ---------------------------------
    function buildPanelE(g)
        addSlider(g, 1, 'Mirror-screen dist.', [0.3 6], P.D, @(v) sprintf('%.2f m', v), ...
            @(v) setParam('D', v, 'geo'));
        addSlider(g, 2, 'Laser brightness', [0.1 1], P.brightness, @(v) sprintf('%.0f %%', 100 * v), ...
            @(v) setParam('brightness', v, 'none'));
        addSlider(g, 3, 'Laser spot size', [1 8], P.spot, @(v) sprintf('%.1f px', v), ...
            @(v) setParam('spot', v, 'none'));
        addSlider(g, 4, 'Projection scale', [log10(0.05) log10(100)], log10(P.viewHalf_cm), ...
            @(v) sprintf('%s%.3g cm', char(177), 10^v), @(v) setParam('viewHalf_cm', 10^v, 'none'));
        addSlider(g, 5, 'Exposure / persistence', [log10(1) log10(200)], log10(P.exposure * 1e3), ...
            @(v) sprintf('%.1f ms', 10^v), @(v) setParam('exposure', 10^v * 1e-3, 'none'));
        addSlider(g, 6, 'Disturbance (base)', [0 0.5], P.noise, @(v) sprintf('%.3f m/s%s', v, char(178)), ...
            @(v) setParam('noise', v, 'none'));
        H.autoFit = uicheckbox(g, 'Text', 'Auto-fit screen view to pattern (sizes still reported in cm)', ...
            'Value', P.autoFit, 'FontColor', CLR.text, ...
            'ValueChangedFcn', @(src, ~) setParam('autoFit', src.Value, 'none'));
        H.autoFit.Layout.Row = 7; H.autoFit.Layout.Column = [1 3];
    end

    % ---------------- Panel F: playback & sweep -----------------------
    function buildPanelF(g)
        bg = uigridlayout(g, [1 5]);
        bg.Layout.Row = [1 2]; bg.Layout.Column = [1 3];
        bg.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};
        bg.Padding = [0 0 0 0]; bg.ColumnSpacing = 4; bg.BackgroundColor = CLR.panel;
        H.playBtn = addButton(bg, 1, 1, [char(9654) ' Play'], @() startPlayback());
        H.playBtn.BackgroundColor = [0.15 0.45 0.2];
        addButton(bg, 1, 2, 'Pause', @() stopPlayback());
        addButton(bg, 1, 3, 'Reset', @() onReset());
        addButton(bg, 1, 4, 'Kick', @() onKick());
        addButton(bg, 1, 5, 'Verify', @() showVerification());
        sg = uigridlayout(g, [1 3]);
        sg.Layout.Row = 3; sg.Layout.Column = [1 3];
        sg.ColumnWidth = {'1x', '1x', '1x'};
        sg.Padding = [0 0 0 0]; sg.ColumnSpacing = 4; sg.BackgroundColor = CLR.panel;
        H.sweepUpBtn = addButton(sg, 1, 1, ['Sweep ' char(9650)], @() startSweep(1));
        H.sweepDnBtn = addButton(sg, 1, 2, ['Sweep ' char(9660)], @() startSweep(-1));
        addButton(sg, 1, 3, 'Stop sweep', @() stopSweep());
        addSlider(g, 4, 'Sweep rate', [log10(0.1) log10(12)], log10(P.sweepRate), ...
            @(v) sprintf('%.2g st/s', 10^v), @(v) setParam('sweepRate', 10^v, 'none'));
        addSlider(g, 5, 'Simulation speed', [-3 0], log10(P.speed), ...
            @(v) sprintf('%.3g%s', 10^v, GL.times), @(v) setParam('speed', 10^v, 'none'));
        H.modeDrop = addDropdown(g, 6, 'Display model', ...
            {'Physics model (speaker-spring-mirror)', 'Idealized Lissajous (comparison)'}, ...
            {'Physics', 'Lissajous'}, P.mode, @(v) setParam('mode', v, 'none'));
        H.show3D = uicheckbox(g, 'Text', 'Update 3D schematic', 'Value', P.show3D, ...
            'FontColor', CLR.text, 'ValueChangedFcn', @(src, ~) setParam('show3D', src.Value, 'none'));
        H.show3D.Layout.Row = 7; H.show3D.Layout.Column = [1 2];
        H.timeLbl = addLabel(g, 8, [1 3], '', CLR.value);
    end

    % ---------------- Panel G: idealized Lissajous --------------------
    function buildPanelG(g)
        addSlider(g, 1, 'x frequency ratio a', [1 9], P.Lis.a, @(v) sprintf('%d', round(v)), ...
            @(v) setParam('Lis.a', round(v), 'none'));
        addSlider(g, 2, 'y frequency ratio b', [1 9], P.Lis.b, @(v) sprintf('%d', round(v)), ...
            @(v) setParam('Lis.b', round(v), 'none'));
        addSlider(g, 3, ['Phase ' GL.delta], [0 360], P.Lis.delta_deg, @(v) sprintf('%.0f%s', v, GL.deg), ...
            @(v) setParam('Lis.delta_deg', v, 'none'));
        addSlider(g, 4, 'Amplitude A (x)', [0.1 30], P.Lis.A_cm, @(v) sprintf('%.1f cm', v), ...
            @(v) setParam('Lis.A_cm', v, 'none'));
        addSlider(g, 5, 'Amplitude B (y)', [0.1 30], P.Lis.B_cm, @(v) sprintf('%.1f cm', v), ...
            @(v) setParam('Lis.B_cm', v, 'none'));
        addSlider(g, 6, 'x detune', [0 5], P.Lis.detune, @(v) sprintf('%.2f Hz', v), ...
            @(v) setParam('Lis.detune', v, 'none'));
        H.bassTreble = uicheckbox(g, 'Text', 'HYPOTHETICAL bass(x)/treble(y) two-frequency view', ...
            'Value', P.Lis.bassTreble, 'FontColor', CLR.warn, ...
            'ValueChangedFcn', @(src, ~) setParam('Lis.bassTreble', src.Value, 'none'));
        H.bassTreble.Layout.Row = 7; H.bassTreble.Layout.Column = [1 3];
    end

    % ---------------- Projection axes ---------------------------------
    function buildProjectionAxes(parent)
        ax = uiaxes(parent);
        ax.Layout.Row = 1; ax.Layout.Column = 1;
        ax.Color = CLR.screen;
        ax.XColor = [0.45 0.45 0.5]; ax.YColor = [0.45 0.45 0.5];
        ax.GridColor = [0.35 0.35 0.4]; ax.GridAlpha = 0.18;
        ax.XGrid = 'on'; ax.YGrid = 'on'; ax.Box = 'on';
        ax.FontSize = 10;
        ax.Title.String = 'Projection screen (laser trajectory, cm)';
        ax.Title.Color = CLR.text;
        ax.XLabel.String = 'screen horizontal (cm)';
        ax.YLabel.String = 'screen vertical (cm)';
        ax.DataAspectRatio = [1 1 1];
        ax.XLim = [-1 1] * P.viewHalf_cm; ax.YLim = [-1 1] * P.viewHalf_cm;
        ax.NextPlot = 'add';
        try, disableDefaultInteractivity(ax); catch, end
        try, ax.Toolbar.Visible = 'off'; catch, end
        ax.BackgroundColor = CLR.bg;
        H.ax = ax;
        % Layered rendering of ONE trajectory: glow -> fading trail -> core -> head
        H.glow = line(ax, NaN, NaN, 'Color', [0.45 0 0], 'LineWidth', 10);
        H.trail = patch(ax, 'Vertices', [NaN NaN], 'Faces', 1, 'FaceVertexCData', [0 0 0], ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
        H.core = line(ax, NaN, NaN, 'Color', [1 0.6 0.55], 'LineWidth', 1);
        H.head = line(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerSize', 8, 'MarkerFaceColor', [1 0.85 0.8], 'MarkerEdgeColor', [1 0.2 0.15]);
        H.info = text(ax, 0.01, 0.985, '', 'Units', 'normalized', 'Color', [0.85 0.85 0.9], ...
            'FontSize', 11, 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'Interpreter', 'none');
        H.warnTxt = text(ax, 0.01, 0.02, '', 'Units', 'normalized', 'Color', CLR.warn, ...
            'FontSize', 11, 'VerticalAlignment', 'bottom', 'FontName', 'Consolas', 'Interpreter', 'none');
    end

    % ---------------- Time-signal tab ---------------------------------
    function buildSignalsTab(tab)
        g = uigridlayout(tab, [3 1]);
        g.BackgroundColor = CLR.bg; g.Padding = [4 4 4 4]; g.RowSpacing = 4;
        H.axSpk = darkAxes(g, 1, 'Speaker cone displacement', 'time (ms)', ['s (' GL.mu 'm)']);
        H.lSpkOut = line(H.axSpk, NaN, NaN, 'Color', [1 0.6 0.2], 'LineWidth', 1.5);
        H.lSpkIn  = line(H.axSpk, NaN, NaN, 'Color', [0.3 0.75 1], 'LineWidth', 1.5);
        legend(H.axSpk, {'outer woofer s_{out}', 'inner tweeter s_{in}'}, 'TextColor', CLR.text, ...
            'Color', CLR.panel, 'EdgeColor', [0.3 0.3 0.3], 'Location', 'northeast');
        H.axSpr = darkAxes(g, 2, 'Spring tip displacement relative to cone', 'time (ms)', ['z (' GL.mu 'm)']);
        H.lSpr = line(H.axSpr, NaN, NaN, 'Color', [0.6 1 0.5], 'LineWidth', 1.5);
        H.axTilt = darkAxes(g, 3, 'Mirror tilt', 'time (ms)', [GL.theta ' (mrad)']);
        H.lThX = line(H.axTilt, NaN, NaN, 'Color', [1 0.35 0.35], 'LineWidth', 1.5);
        H.lThY = line(H.axTilt, NaN, NaN, 'Color', [1 0.85 0.3], 'LineWidth', 1.5);
        legend(H.axTilt, {[GL.theta '_x horizontal tilt (about vertical axis)'], ...
            [GL.theta '_y vertical tilt (about horizontal axis)']}, 'TextColor', CLR.text, ...
            'Color', CLR.panel, 'EdgeColor', [0.3 0.3 0.3], 'Location', 'northeast');
    end

    function ax = darkAxes(g, row, ttl, xl, yl)
        ax = uiaxes(g);
        ax.Layout.Row = row; ax.Layout.Column = 1;
        ax.Color = [0.06 0.065 0.08];
        ax.XColor = CLR.text; ax.YColor = CLR.text; ax.ZColor = CLR.text;
        ax.GridColor = [0.5 0.5 0.55]; ax.GridAlpha = 0.2;
        ax.XGrid = 'on'; ax.YGrid = 'on'; ax.Box = 'on';
        ax.FontSize = 9;
        ax.Title.String = ttl; ax.Title.Color = CLR.text;
        ax.XLabel.String = xl; ax.YLabel.String = yl;
        ax.BackgroundColor = CLR.bg;
        ax.NextPlot = 'add';
        try, disableDefaultInteractivity(ax); catch, end
    end

    % ---------------- Resonance tab -----------------------------------
    function buildResonanceTab(tab)
        g = uigridlayout(tab, [2 1]);
        g.RowHeight = {'1x', 76};
        g.BackgroundColor = CLR.bg; g.Padding = [4 4 4 4];
        H.axRes = darkAxes(g, 1, 'Small-signal (linear) response to a sine at the current amplifier settings', ...
            'drive frequency (Hz)', 'amplitude');
        H.axRes.XScale = 'log'; H.axRes.YScale = 'log';
        H.axRes.XLim = [1 P.fMax];
        H.lResTx = line(H.axRes, NaN, NaN, 'Color', [1 0.35 0.35], 'LineWidth', 1.5);
        H.lResTy = line(H.axRes, NaN, NaN, 'Color', [1 0.85 0.3], 'LineWidth', 1.5);
        H.lResZ  = line(H.axRes, NaN, NaN, 'Color', [0.6 1 0.5], 'LineWidth', 1.2);
        H.lResSpot = line(H.axRes, NaN, NaN, 'Color', [0.8 0.5 1], 'LineWidth', 1.2, 'LineStyle', '--');
        H.lResModes = line(H.axRes, NaN, NaN, 'Color', [0.6 0.6 0.7], 'LineStyle', ':', 'LineWidth', 1);
        H.lResNow = line(H.axRes, NaN, NaN, 'Color', [1 1 1], 'LineWidth', 1.2);
        legend(H.axRes, {[GL.theta '_x (mrad)'], [GL.theta '_y (mrad)'], ['z (' GL.mu 'm)'], ...
            'spot half-size (cm)', 'coupled natural freqs', 'current f'}, ...
            'TextColor', CLR.text, 'Color', CLR.panel, 'EdgeColor', [0.3 0.3 0.3], 'Location', 'southwest');
        H.resNote = uilabel(g, 'Text', '', 'FontColor', CLR.text, 'WordWrap', 'on', 'FontSize', 11);
        H.resNote.Layout.Row = 2; H.resNote.Layout.Column = 1;
    end

    % ---------------- 3D schematic tab --------------------------------
    function buildSchematicTab(tab)
        g = uigridlayout(tab, [2 4]);
        g.RowHeight = {24, '1x'};
        g.ColumnWidth = {50, '1x', 90, '1x'};
        g.BackgroundColor = CLR.bg; g.Padding = [4 4 4 4];
        lb = uilabel(g, 'Text', 'View', 'FontColor', CLR.text); lb.Layout.Row = 1; lb.Layout.Column = 1;
        H.viewDrop = uidropdown(g, 'Items', {'Apparatus close-up', 'Full ray path'}, ...
            'BackgroundColor', CLR.field, 'FontColor', CLR.text, ...
            'ValueChangedFcn', @(~, ~) applySchematicView());
        H.viewDrop.Layout.Row = 1; H.viewDrop.Layout.Column = 2;
        lb = uilabel(g, 'Text', 'Tilt drawn', 'FontColor', CLR.text); lb.Layout.Row = 1; lb.Layout.Column = 3;
        H.exagDrop = uidropdown(g, 'Items', {[GL.times '1 (true)'], [GL.times '10'], [GL.times '30']}, ...
            'ItemsData', [1 10 30], 'Value', 1, ...
            'BackgroundColor', CLR.field, 'FontColor', CLR.text);
        H.exagDrop.Layout.Row = 1; H.exagDrop.Layout.Column = 4;
        ax = uiaxes(g);
        ax.Layout.Row = 2; ax.Layout.Column = [1 4];
        ax.Color = [0.07 0.075 0.09];
        ax.XColor = CLR.text; ax.YColor = CLR.text; ax.ZColor = CLR.text;
        ax.GridColor = [0.5 0.5 0.55]; ax.GridAlpha = 0.15;
        ax.XGrid = 'on'; ax.YGrid = 'on'; ax.ZGrid = 'on';
        ax.BackgroundColor = CLR.bg;
        ax.XLabel.String = 'x (m)'; ax.YLabel.String = 'y (m, speaker axis)'; ax.ZLabel.String = 'z (m, up)';
        ax.Title.String = 'Apparatus schematic (same geometry as the ray tracer)';
        ax.Title.Color = CLR.text;
        ax.DataAspectRatio = [1 1 1];
        ax.NextPlot = 'add';
        view(ax, [-38 18]);
        H.ax3 = ax;
    end

    % -----------------------------------------------------------------
    % Parameter handling
    % -----------------------------------------------------------------
    function setParam(name, value, kind)
        parts = strsplit(name, '.');
        if numel(parts) == 1
            P.(name) = value;
        else
            P.(parts{1}).(parts{2}) = value;
        end
        switch kind
            case 'dyn'      % changes the mechanical / signal model
                S.needRebuild = true;
                S.resDirty = true;
            case 'geo'      % changes only the optical geometry
                S.geoDirty = true;
                S.schemDirty = true;
        end
        if strcmp(name, 'autoFit') && ~value
            S.viewHalf = P.viewHalf_cm;
        end
        if any(strcmp(name, {'e_h', 'e_v'}))
            S.schemDirty = true;    % spring/mirror drawing depends on the offsets
        end
        if ~S.playing
            processPending();
            renderFrame(true);
        end
    end

    function processPending()
        % Model rebuilds keep the mechanical state and the tone phase, so
        % any change is followed by the physical (integrated) response.
        if S.needRebuild
            St = buildStructure(P);
            Mdl = buildModel(P, St);
            G = buildGeometry(P);
            S.needRebuild = false;
            S.resDirty = true;
            updateReadouts();
        end
        if S.geoDirty
            G = buildGeometry(P);
            S.geoDirty = false;
            S.resDirty = true;
        end
    end

    function resetSimulation()
        % Tone running at t = 0 with the mirror on its small-signal steady
        % state; the nonlinear integration then takes over.
        St = buildStructure(P);
        Mdl = buildModel(P, St);
        G = buildGeometry(P);
        S.t = 0;
        S.phase = 0;
        S.an = 0;
        S.simDebt = 0;
        S.x = ssStatePhase(Mdl, 0);
        S.ampEst = max(abs(S.x(2:3)));
        S.ampEst = max(S.ampEst, sum(max(abs(Mdl.QF(2:3, :)), [], 1)));
        % pre-fill the persistence buffer with the same steady state
        S.buf(:) = NaN; S.bufHead = 0; S.bufCount = 0;
        dt0 = stepSize(Mdl, P, S.ampEst);
        nPre = min(floor(P.exposure / dt0) + 1, floor(S.Nb / 2));
        tPre = -(nPre - 1:-1:0) * dt0;
        phPre = mod(Mdl.f * tPre, 1);
        bufAppend([tPre; tableEval(Mdl.qL, phPre); tableEval(Mdl.sL, phPre)]);
        S.needRebuild = false;
        S.geoDirty = false;
        S.resDirty = true;
        S.schemDirty = true;
    end

    function onReset()
        stopSweep();
        resetSimulation();
        updateReadouts();
        renderFrame(true);
    end

    function onKick()
        % A tap on the speaker: an angular velocity impulse. The free
        % ringing beats with the forced motion and, with contact hardening,
        % elliptical orbits precess, so the figure visibly turns.
        processPending();
        wk = 2 * pi * St.fn(min(2, numel(St.fn)));
        S.x(5) = S.x(5) + P.theta0 * wk;
        S.x(6) = S.x(6) - 0.7 * P.theta0 * wk;
        S.ampEst = max(S.ampEst, 1.5 * P.theta0);
        if ~S.playing, startPlayback(); end
    end

    % -----------------------------------------------------------------
    % Frequency controls
    % -----------------------------------------------------------------
    function setFrequency(f, source)
        if ~isfinite(f), return; end
        if ~strcmp(source, 'sweep'), stopSweep(); end
        f = round(f * 1000) / 1000;
        f = min(max(f, P.fMin), P.fMax);
        P.f = f;
        S.needRebuild = true;
        S.resDirty = true;
        if strcmp(source, 'sweep')
            S.freqWidgetsDirty = true;
        else
            refreshFrequencyWidgets(source);
        end
        if ~S.playing
            processPending();
            renderFrame(true);
        end
    end

    function refreshFrequencyWidgets(source)
        f = P.f;
        if ~strcmp(source, 'slider'), H.freqSlider.Value = log10(f); end
        if ~strcmp(source, 'edit'), H.freqEdit.Value = f; end
        H.freqValLbl.Text = sprintf('%s Hz', formatFreq(f));
        [nm, cents, midi] = nearestNote(f);
        H.noteInfo.Text = sprintf('Nearest note: %s  (%+.1f cents)   f(n) = 440*2^((n-69)/12)', nm, cents);
        if ~strcmp(source, 'note')
            if abs(cents) < 0.05 && midi >= 21 && midi <= 119
                H.noteDrop.Value = midi;
            else
                H.noteDrop.Value = 0;
            end
        end
        H.freqBig.Text = sprintf('%s Hz', formatFreq(f));
    end

    function onNoteSelected(m)
        if m > 0
            setFrequency(noteFrequency(m), 'note');
        end
    end

    function setStep(v)
        P.fStep = v;
    end

    function startSweep(direction)
        S.sweepDir = direction;
        if direction > 0
            H.sweepUpBtn.BackgroundColor = [0.15 0.45 0.2];
            H.sweepDnBtn.BackgroundColor = CLR.button;
        else
            H.sweepDnBtn.BackgroundColor = [0.15 0.45 0.2];
            H.sweepUpBtn.BackgroundColor = CLR.button;
        end
        if ~S.playing, startPlayback(); end
    end

    function stopSweep()
        S.sweepDir = 0;
        if isfield(H, 'sweepUpBtn') && isvalid(H.sweepUpBtn)
            H.sweepUpBtn.BackgroundColor = CLR.button;
            H.sweepDnBtn.BackgroundColor = CLR.button;
        end
    end

    % -----------------------------------------------------------------
    % Playback
    % -----------------------------------------------------------------
    function startPlayback()
        if S.playing || isempty(S.timer) || ~isvalid(S.timer), return; end
        S.playing = true;
        S.lastTic = tic;
        H.playBtn.Text = [char(9654) ' Playing'];
        start(S.timer);
    end

    function stopPlayback()
        if ~isempty(S.timer) && isvalid(S.timer) && strcmp(S.timer.Running, 'on')
            stop(S.timer);
        end
        S.playing = false;
        if isfield(H, 'playBtn') && isvalid(H.playBtn)
            H.playBtn.Text = [char(9654) ' Play'];
        end
    end

    function onTick()
        try
            if ~isvalid(H.fig), return; end
            dtWall = toc(S.lastTic);
            S.lastTic = tic;
            dtWall = min(dtWall, 0.1);          % avoid jumps after GUI stalls
            dtSim = dtWall * P.speed;
            if S.sweepDir ~= 0
                % sweep rate is in semitones per SIMULATED second
                fNew = P.f * 2^(S.sweepDir * P.sweepRate / 12 * dtSim);
                if fNew >= P.fMax || fNew <= P.fMin
                    fNew = min(max(fNew, P.fMin), P.fMax);
                    stopSweep();
                end
                setFrequency(fNew, 'sweep');
            end
            processPending();
            advance(dtSim);
            if S.freqWidgetsDirty && mod(S.frame, 3) == 0
                refreshFrequencyWidgets('sweep');
                S.freqWidgetsDirty = false;
            end
            S.effSpeed = S.lastSteps * S.lastDt / max(dtWall, eps);
            renderFrame(false);
        catch err
            stopPlayback();
            if isvalid(H.fig)
                uialert(H.fig, sprintf('Animation stopped: %s', err.message), 'Simulation error');
            end
        end
    end

    function advance(dtSim)
        % Integrate the nonlinear model over the newly elapsed simulated
        % time only. A step budget keeps the GUI responsive; if it is hit,
        % simulated time runs slower than requested (reported on screen).
        dt = stepSize(Mdl, P, S.ampEst);
        S.simDebt = S.simDebt + dtSim;
        n = floor(S.simDebt / dt);
        S.cpuLimited = n > S.stepBudget;
        if S.cpuLimited
            n = S.stepBudget;
            S.simDebt = 0;
        else
            S.simDebt = S.simDebt - n * dt;
        end
        S.lastSteps = n; S.lastDt = dt;
        if n < 1, return; end
        rho = exp(-2 * pi * P.noiseCorner * dt);
        xi = P.noise * sqrt(1 - rho^2) * randn(S.rng, n, 1);
        tic0 = tic;
        [xNew, phNew, anNew, rec] = integrateChunk(S.x, S.phase, S.an, n, dt, Mdl, rho, xi);
        el = toc(tic0);
        if ~all(isfinite(xNew)) || ~all(isfinite(rec(:)))
            % numerical safeguard: restart from the small-signal steady state
            S.x = ssStatePhase(Mdl, S.phase);
            S.an = 0;
            S.nanReset = true;
            return
        end
        tRec = S.t + (1:n) * dt;
        S.t = tRec(end);
        S.x = xNew; S.phase = phNew; S.an = anNew;
        S.ampEst = max(abs(rec(3:4, :)), [], 'all');
        bufAppend([tRec; rec(2:4, :); tableEval(Mdl.sL, rec(1, :))]);
        % adapt the step budget to keep integration under ~12 ms per frame
        if el > 0
            S.stepBudget = round(min(60000, max(800, 0.7 * S.stepBudget + 0.3 * n * 0.012 / el)));
        end
    end

    function onClose()
        try
            if ~isempty(S.timer) && isvalid(S.timer)
                stop(S.timer);
                delete(S.timer);
            end
        catch
        end
        try
            if ~isempty(S.verifyFig) && isvalid(S.verifyFig), delete(S.verifyFig); end
        catch
        end
        delete(H.fig);
    end

    % -----------------------------------------------------------------
    % Ring buffer of simulated samples (the visible persistence trail)
    % -----------------------------------------------------------------
    function bufAppend(D)
        n = size(D, 2);
        if n == 0, return; end
        if n >= S.Nb
            D = D(:, end - S.Nb + 1:end);
            n = S.Nb;
        end
        idx = mod(S.bufHead + (0:n - 1), S.Nb) + 1;
        S.buf(:, idx) = D;
        S.bufHead = idx(end);
        S.bufCount = min(S.Nb, S.bufCount + n);
    end

    function D = bufRecent(tFrom)
        if S.bufCount == 0
            D = zeros(6, 0);
            return
        end
        idx = mod(S.bufHead - S.bufCount + (0:S.bufCount - 1), S.Nb) + 1;
        tt = S.buf(1, idx);
        first = find(tt >= tFrom, 1);
        if isempty(first)
            D = S.buf(:, idx(end));
        else
            D = S.buf(:, idx(first:end));
        end
    end

    function f = fastestContent()
        % highest frequency present in the motion (drive harmonics or modes)
        f = max(Mdl.f * Mdl.nmaxF, Mdl.fnMax * sqrt(1 + 3 * P.eta_nl * (S.ampEst / P.theta_ref)^2));
    end

    % -----------------------------------------------------------------
    % Rendering
    % -----------------------------------------------------------------
    function renderFrame(forceAll)
        if nargin < 1, forceAll = false; end
        tStart = tic;
        S.frame = S.frame + 1;
        truncated = false;
        tiltWarn = false;

        % ---- the visible trajectory (rolling exposure window) ----
        if strcmp(P.mode, 'Physics')
            D = bufRecent(S.t - P.exposure);
            [D, truncated] = selectSamples(D, S.NsMax, fastestContent());
            [X, Y, tiltWarn] = qToScreen(D(2:4, :), G, P);
            tv = D(1, :);
        else
            [X, Y, tv] = lissajousTrail(S.t, P, S.NsMax);
        end
        Ns = numel(X);
        if Ns == 0
            X = NaN; Y = NaN; tv = S.t; Ns = 1;
        end

        % ---- update the ONE laser trajectory (existing objects only) ----
        br = P.brightness;
        age = linspace(0, 1, Ns)';
        red = br * (0.10 + 0.90 * age .^ 0.8);
        cols = [red, 0.07 * red, 0.06 * red];
        set(H.trail, 'Vertices', [X(:) Y(:); NaN NaN], 'Faces', 1:(Ns + 1), ...
            'FaceVertexCData', [cols; 0 0 0], 'LineWidth', P.spot);
        set(H.glow, 'XData', X, 'YData', Y, 'Color', br * [0.42 0.0 0.02], ...
            'LineWidth', P.spot * 3.2);
        nCore = max(2, round(0.3 * Ns));
        idxCore = max(1, Ns - nCore + 1):Ns;
        set(H.core, 'XData', X(idxCore), 'YData', Y(idxCore), ...
            'Color', min(1, [1 0.55 0.5] * (0.4 + 0.6 * br)), 'LineWidth', max(0.5, 0.35 * P.spot));
        set(H.head, 'XData', X(end), 'YData', Y(end), 'MarkerSize', 2 + 2.2 * P.spot);

        % ---- view limits -------------------------------------------------
        ext = max(abs([X(:); Y(:)]), [], 'omitnan');
        if isempty(ext) || ~isfinite(ext), ext = 0; end
        if P.autoFit
            target = max(1.35 * ext, 0.01);
            if target > S.viewHalf
                S.viewHalf = target;
            elseif target < 0.45 * S.viewHalf
                S.viewHalf = S.viewHalf + 0.25 * (target - S.viewHalf);
            end
        else
            S.viewHalf = P.viewHalf_cm;
        end
        lim = [-1 1] * S.viewHalf;
        if abs(H.ax.XLim(2) - lim(2)) > 1e-3 * lim(2)
            H.ax.XLim = lim; H.ax.YLim = lim;
        end

        % ---- pattern statistics (comparison with photographs) -----------
        st = patternStats(X, Y);
        [nm, cents] = nearestNote(P.f);
        if strcmp(P.mode, 'Physics')
            modeStr = 'PHYSICS MODEL (nonlinear, integrated in time)';
        elseif P.Lis.bassTreble
            modeStr = 'HYPOTHETICAL two-frequency (NOT how bass/treble work)';
        else
            modeStr = 'IDEALIZED LISSAJOUS (not the physical system)';
        end
        sweepStr = '';
        if S.sweepDir > 0, sweepStr = sprintf('  SWEEP UP %.2g st/s', P.sweepRate); end
        if S.sweepDir < 0, sweepStr = sprintf('  SWEEP DOWN %.2g st/s', P.sweepRate); end
        H.info.String = sprintf(['%s Hz  (%s %+.0f c)  %s wave%s\n%s\n' ...
            'pattern %.2f x %.2f cm, major axis %+.0f deg   view %s%.3g cm'], ...
            formatFreq(P.f), nm, cents, P.wave, sweepStr, modeStr, st.len, st.wid, st.ang, ...
            char(177), S.viewHalf);
        warn = {};
        if S.cpuLimited && S.playing
            warn{end + 1} = sprintf('CPU-limited: simulated time running at %.2g%s instead of %.2g%s', ...
                S.effSpeed, GL.times, P.speed, GL.times);
        end
        if truncated && strcmp(P.mode, 'Physics')
            warn{end + 1} = sprintf('High frequency: last %.3g ms of the exposure drawn (anti-aliased)', ...
                1e3 * (tv(end) - tv(1)));
        end
        if tiltWarn
            warn{end + 1} = 'Tilt beyond 0.25 rad: outside the small-tilt contact model (mirror would rattle)';
        end
        if S.nanReset
            warn{end + 1} = 'Numerical safeguard: state restarted from steady state';
            S.nanReset = false;
        end
        if ~S.playing && strcmp(P.mode, 'Physics')
            warn{end + 1} = 'Paused: the trail shows the motion before any change. Press Play.';
        end
        H.warnTxt.String = strjoin(warn, newline);
        H.timeLbl.Text = sprintf('t_sim = %.4f s  steps/frame %d  dt %.3g %ss', ...
            S.t, S.lastSteps, 1e6 * S.lastDt, GL.mu);
        H.status.Text = strsplit(statusText(st), newline);

        % ---- secondary plots (throttled) ---------------------------------
        if forceAll || mod(S.frame, 3) == 0
            updateSignalPlots();
        end
        if P.show3D && (forceAll || mod(S.frame, 4) == 0)
            if S.schemDirty, drawSchematicStatic(); end
            updateSchematic();
        end
        if S.resDirty && (forceAll || mod(S.frame, 10) == 0)
            updateResonancePlot();
        end
        drawnow limitrate;

        % ---- performance safeguard: adapt drawn-sample budget ----------
        el = toc(tStart);
        if el > 0.020
            S.NsMax = max(600, round(S.NsMax * 0.85));
        elseif el < 0.008
            S.NsMax = min(3000, round(S.NsMax * 1.05) + 1);
        end
    end

    function txt = statusText(st)
        Q1 = Mdl.QF(:, 1);
        dphi = angle(Q1(2) / (Q1(3) + eps)) * 180 / pi;
        txt = sprintf([' Small-signal fundamental: |%s_x| = %.3g mrad, |%s_y| = %.3g mrad, phase(%s_y - %s_x) = %+.0f%s' ...
            '\n Trail: %.2f x %.2f cm (axis ratio %.2f)   max tilt %.3g mrad   harmonics in drive: %d'], ...
            GL.theta, 1e3 * abs(Q1(3)), GL.theta, 1e3 * abs(Q1(2)), GL.theta, GL.theta, dphi, GL.deg, ...
            st.len, st.wid, st.wid / max(st.len, eps), 1e3 * S.ampEst, numel(Mdl.nF));
    end

    function updateSignalPlots()
        Wp = min(max(3 / P.f, 0.02), 0.5);
        D = bufRecent(S.t - Wp);
        D = selectSamples(D, 2000, fastestContent());
        if size(D, 2) < 2, return; end
        tm = (D(1, :) - S.t) * 1e3;
        set(H.lSpkOut, 'XData', tm, 'YData', D(5, :) * 1e6);
        set(H.lSpkIn,  'XData', tm, 'YData', D(6, :) * 1e6);
        set(H.lSpr,    'XData', tm, 'YData', D(2, :) * 1e6);
        set(H.lThX,    'XData', tm, 'YData', D(4, :) * 1e3);
        set(H.lThY,    'XData', tm, 'YData', D(3, :) * 1e3);
        xl = [tm(1) max(tm(end), tm(1) + 1e-6)];
        H.axSpk.XLim = xl; H.axSpr.XLim = xl; H.axTilt.XLim = xl;
    end

    function updateResonancePlot()
        S.resDirty = false;
        fs = logspace(0, log10(P.fMax), 500);
        U = (P.volume / 100) * eqResponse(fs, P) .* (-1i);
        Q = harmonicResponse(St, P, fs, U);
        floorV = 1e-9;
        thx = max(abs(Q(3, :)) * 1e3, floorV);
        thy = max(abs(Q(2, :)) * 1e3, floorV);
        zz  = max(abs(Q(1, :)) * 1e6, floorV);
        spot = max(2 * P.D * hypot(abs(Q(2, :)), abs(Q(3, :))) * 100, floorV);
        set(H.lResTx, 'XData', fs, 'YData', thx);
        set(H.lResTy, 'XData', fs, 'YData', thy);
        set(H.lResZ,  'XData', fs, 'YData', zz);
        set(H.lResSpot, 'XData', fs, 'YData', spot);
        allY = [thx thy zz spot];
        yTop = max(allY) * 3;
        yBot = max(min(allY), yTop * 1e-9);
        fm = St.fn(:)';
        xs = reshape([fm; fm; nan(size(fm))], 1, []);
        ys = reshape([yBot * ones(size(fm)); yTop * ones(size(fm)); nan(size(fm))], 1, []);
        set(H.lResModes, 'XData', xs, 'YData', ys);
        set(H.lResNow, 'XData', [P.f P.f], 'YData', [yBot yTop]);
        H.axRes.YLim = [yBot yTop];
        fmL = sprintf('%.1f, ', St.fn); fmL = fmL(1:end - 2);
        fr = St.fn(:)';
        sub = sprintf('%.0f, ', sort([fr / 2, fr / 3])); sub = sub(1:end - 2);
        H.resNote.Text = sprintf([' Coupled natural frequencies: %s Hz (undamped). Lines and ellipses ' ...
            'turn as f passes each mode. With harmonic distortion, multi-lobe figures appear near f_n/2 ' ...
            'and f_n/3 = %s Hz. Contact hardening bends the peaks to the right at large amplitude, so ' ...
            'up- and down-sweeps can jump at different frequencies.'], fmL, sub);
    end

    function updateReadouts()
        if ~isfield(H, 'eqInfo'), return; end
        Heq = eqResponse(P.f, P);
        H.eqInfo.Text = sprintf('EQ gain at %s Hz: %+.2f dB   (linear volume %.2f)', ...
            formatFreq(P.f), 20 * log10(abs(Heq)), P.volume / 100);
        H.spkInfo.Text = sprintf('At f: |s_out| = %.3g %sm, |s_in| = %.3g %sm (fundamental)', ...
            1e6 * abs(Mdl.So1), GL.mu, 1e6 * abs(Mdl.Si1), GL.mu);
        H.fnInfo1.Text = sprintf('1-DOF f_n = %.1f Hz (k_eff = %.0f N/m, m = %.2f g)', ...
            St.fn1, St.keff1, P.m * 1e3);
        H.fnInfo2.Text = sprintf('Coupled modes: %.1f, %.1f, %.1f Hz (tilt pair split %.0f%%)', ...
            St.fn(1), St.fn(2), St.fn(3), 100 * (St.fn(3) / St.fn(2) - 1));
        [dmin, im] = min(abs(P.f - St.fn) ./ St.fn);
        if dmin < 2 * P.zeta
            H.resLamp.Color = [1 0.2 0.2];
            H.resText.Text = sprintf('AT RESONANCE with mode %d (%.1f Hz)', im, St.fn(im));
        elseif dmin < 10 * P.zeta
            H.resLamp.Color = [1 0.8 0.2];
            H.resText.Text = sprintf('Near mode %d (%.1f Hz), %.0f%% away', im, St.fn(im), 100 * dmin);
        else
            H.resLamp.Color = [0.2 0.8 0.3];
            H.resText.Text = sprintf('Off resonance (nearest mode %.1f Hz)', St.fn(im));
        end
    end

    % -----------------------------------------------------------------
    % 3D schematic
    % -----------------------------------------------------------------
    function drawSchematicStatic()
        S.schemDirty = false;
        ax = H.ax3;
        cla(ax);
        ax.NextPlot = 'add';
        wood = [0.80 0.63 0.40];
        zc = P.geo.speakerCenter(3);
        % wooden base (horizontal) and speaker panel (vertical)
        patch(ax, 'XData', [-0.16 0.16 0.16 -0.16], 'YData', [-0.03 -0.03 0.32 0.32], ...
            'ZData', [0 0 0 0], 'FaceColor', wood, 'EdgeColor', [0.5 0.4 0.25]);
        patch(ax, 'XData', [-0.16 0.16 0.16 -0.16], 'YData', [-0.03 -0.03 -0.03 -0.03], ...
            'ZData', [0 0 0.34 0.34], 'FaceColor', wood * 0.9, 'EdgeColor', [0.5 0.4 0.25]);
        % woofer surround and cone
        th = linspace(0, 2 * pi, 60);
        [TT, RR] = meshgrid(th, [0.075 0.088]);
        surf(ax, RR .* cos(TT), zeros(size(TT)), zc + RR .* sin(TT), ...
            'FaceColor', [0.12 0.12 0.12], 'EdgeColor', 'none');
        [TT, RR] = meshgrid(th, linspace(0.022, 0.075, 8));
        surf(ax, RR .* cos(TT), coneY(RR), zc + RR .* sin(TT), ...
            'FaceColor', [0.85 0.47 0.22], 'EdgeColor', 'none', 'FaceAlpha', 0.95);
        % tweeter (inner speaker)
        [TT, YY] = meshgrid(th, [coneY(0.022) 0.006]);
        surf(ax, 0.022 * cos(TT), YY, zc + 0.022 * sin(TT), 'FaceColor', [0.15 0.15 0.15], 'EdgeColor', 'none');
        patch(ax, 'XData', 0.022 * cos(th), 'YData', 0.006 * ones(size(th)), ...
            'ZData', zc + 0.022 * sin(th), 'FaceColor', [0.25 0.25 0.27], 'EdgeColor', 'none');
        % tape, spring (1 inch), mirror pivot
        M0 = G.M0;
        pivot = M0 - [P.e_h; 0; P.e_v];
        root = pivot - [P.L; 0; 0];
        ys = M0(2) - 0.002;
        tx = root(1) + [-0.008 0.006 0.006 -0.008];
        tz = root(3) + [-0.012 -0.012 0.012 0.012];
        patch(ax, 'XData', tx, 'YData', coneY(abs(root(1))) * ones(1, 4) + 0.0005, 'ZData', tz, ...
            'FaceColor', [0.2 0.45 0.95], 'EdgeColor', 'none');
        nTurn = 14;
        sp = linspace(0, 1, 400);
        plot3(ax, root(1) + P.L * sp, ys + 0.0013 * cos(2 * pi * nTurn * sp), ...
            root(3) + 0.0013 * sin(2 * pi * nTurn * sp), 'Color', [0.8 0.8 0.85], 'LineWidth', 1.2);
        % laser pointer on its mount
        L0 = G.L0;
        body = [L0 - 0.07 * G.din, L0];
        plot3(ax, body(1, :), body(2, :), body(3, :), 'Color', [0.2 0.2 0.22], 'LineWidth', 7);
        plot3(ax, [L0(1) L0(1)], [L0(2) L0(2)], [0 L0(3)], 'Color', [0.35 0.35 0.35], 'LineWidth', 3);
        % projection screen (plane perpendicular to the nominal reflected ray)
        hw = 0.45; hh = 0.32;
        cs = [G.S0 + hw * G.u + hh * G.v, G.S0 - hw * G.u + hh * G.v, ...
              G.S0 - hw * G.u - hh * G.v, G.S0 + hw * G.u - hh * G.v];
        patch(ax, 'XData', cs(1, :), 'YData', cs(2, :), 'ZData', cs(3, :), ...
            'FaceColor', [0.78 0.78 0.8], 'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.85);
        text(ax, G.S0(1), G.S0(2), G.S0(3) + hh + 0.04, 'projection screen', 'Color', CLR.text, ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
        text(ax, L0(1), L0(2), L0(3) + 0.03, 'laser', 'Color', CLR.text, 'FontSize', 9);
        text(ax, root(1) - 0.01, M0(2) + 0.01, root(3) - 0.025, 'tape + 1" spring', ...
            'Color', CLR.text, 'FontSize', 9);
        % dynamic objects
        H.s3Mirror = patch(ax, 'XData', NaN, 'YData', NaN, 'ZData', NaN, ...
            'FaceColor', [0.85 0.9 1], 'EdgeColor', [0.4 0.4 0.5], 'FaceAlpha', 0.95);
        H.s3In  = plot3(ax, NaN, NaN, NaN, 'Color', [1 0.1 0.1], 'LineWidth', 1.5);
        H.s3Out = plot3(ax, NaN, NaN, NaN, 'Color', [1 0.25 0.2], 'LineWidth', 1.5);
        H.s3Spot = plot3(ax, NaN, NaN, NaN, 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', 'none');
        H.s3Normal = plot3(ax, NaN, NaN, NaN, 'Color', [0.4 1 0.4], 'LineWidth', 1);
        applySchematicView();
    end

    function applySchematicView()
        ax = H.ax3;
        if strcmp(H.viewDrop.Value, 'Apparatus close-up')
            ax.XLim = [-0.17 0.17]; ax.YLim = [-0.04 0.33]; ax.ZLim = [0 0.35];
        else
            pts = [G.L0, G.M0, G.S0 + 0.5 * G.u + 0.4 * G.v, G.S0 - 0.5 * G.u - 0.4 * G.v, ...
                [-0.17; -0.04; 0], [0.17; 0.33; 0.35]];
            lo = min(pts, [], 2) - 0.05; hi = max(pts, [], 2) + 0.05;
            ax.XLim = [lo(1) hi(1)]; ax.YLim = [lo(2) hi(2)]; ax.ZLim = [min(lo(3), 0) hi(3)];
        end
    end

    function updateSchematic()
        if ~isfield(H, 's3Mirror') || ~isvalid(H.s3Mirror), return; end
        exag = H.exagDrop.Value;
        if strcmp(P.mode, 'Physics')
            q = S.x(1:3);
        else
            q = zeros(3, 1);
        end
        thV = min(max(q(2), -P.tiltMax), P.tiltMax) * exag;
        thH = min(max(q(3), -P.tiltMax), P.tiltMax) * exag;
        [~, ~, valid, hit, dout, n, Pscr] = laserSpot(G, q(1), thV, thH);
        % mirror disc with the current (optionally exaggerated) normal
        [e1, e2] = planeBasis(n);
        ang = linspace(0, 2 * pi, 40);
        c = G.M0 + [0; q(1) + P.e_h * thH - P.e_v * thV; 0];
        disc = c + P.R_mirror * (e1 * cos(ang) + e2 * sin(ang));
        set(H.s3Mirror, 'XData', disc(1, :), 'YData', disc(2, :), 'ZData', disc(3, :));
        set(H.s3In, 'XData', [G.L0(1) hit(1)], 'YData', [G.L0(2) hit(2)], 'ZData', [G.L0(3) hit(3)]);
        nn = c + 0.03 * n;
        set(H.s3Normal, 'XData', [c(1) nn(1)], 'YData', [c(2) nn(2)], 'ZData', [c(3) nn(3)]);
        if valid
            set(H.s3Out, 'XData', [hit(1) Pscr(1)], 'YData', [hit(2) Pscr(2)], 'ZData', [hit(3) Pscr(3)]);
            set(H.s3Spot, 'XData', Pscr(1), 'YData', Pscr(2), 'ZData', Pscr(3));
        else
            e = hit + 0.3 * dout;
            set(H.s3Out, 'XData', [hit(1) e(1)], 'YData', [hit(2) e(2)], 'ZData', [hit(3) e(3)]);
            set(H.s3Spot, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
        end
        if exag > 1
            H.ax3.Title.String = sprintf('Schematic: tilt drawn %s%d (projection uses true angles)', GL.times, exag);
        else
            H.ax3.Title.String = 'Apparatus schematic (same geometry as the ray tracer)';
        end
    end

    % -----------------------------------------------------------------
    % Verification
    % -----------------------------------------------------------------
    function showVerification()
        wasPlaying = S.playing;
        stopPlayback();
        d = uiprogressdlg(H.fig, 'Title', 'Verification', 'Message', 'Running numerical checks...', ...
            'Indeterminate', 'on');
        R = runVerification(P);
        close(d);
        printVerification(R);
        if isempty(S.verifyFig) || ~isvalid(S.verifyFig)
            S.verifyFig = uifigure('Name', 'Numerical verification report', ...
                'Position', [120 120 1180 560], 'Color', CLR.bg);
            g = uigridlayout(S.verifyFig, [2 1]);
            g.RowHeight = {40, '1x'}; g.BackgroundColor = CLR.bg;
            H.vSummary = uilabel(g, 'Text', '', 'FontColor', CLR.text, 'FontSize', 14, 'FontWeight', 'bold');
            H.vTable = uitable(g, 'ColumnName', {'#', 'Test', 'Result', 'Numerical error / metric', 'Details'}, ...
                'ColumnWidth', {30, 290, 60, 160, 'auto'}, 'RowName', {});
        end
        data = cell(numel(R), 5);
        for i = 1:numel(R)
            data(i, :) = {R(i).id, R(i).name, passFail(R(i).pass), sprintf('%.3g', R(i).metric), R(i).detail};
        end
        H.vTable.Data = data;
        try
            for i = 1:numel(R)
                if R(i).pass, colr = [0.75 0.95 0.75]; else, colr = [1 0.75 0.75]; end
                addStyle(H.vTable, uistyle('BackgroundColor', colr), 'cell', [i 3]);
            end
        catch
        end
        H.vSummary.Text = sprintf(['%d / %d checks passed. These establish internal numerical ' ...
            'consistency only; they do not prove agreement with the physical experiment.'], ...
            sum([R.pass]), numel(R));
        figure(S.verifyFig);
        if wasPlaying, startPlayback(); end
    end

end   % ===== end of main function (nested functions above) =====

% =====================================================================
%                     LOCAL FUNCTIONS (PHYSICS CORE)
% =====================================================================

function P = defaultParams()
% All values are illustrative ESTIMATES chosen to give a recognisable
% pattern; replace them with measured values when available.
% ---- tone generator ----
P.f = 440;                 % Hz (A4, as on the tone generator start page)
P.wave = 'Sine';
P.fStep = 1;               % Hz per arrow-button press
P.fMin = 1; P.fMax = 20154;
P.audioBandwidth = 22050;  % Hz: harmonics above this are not generated (44.1 kHz audio)
P.nHarmMax = 512;          % max Fourier harmonics (reaches the mirror modes even for a 1 Hz square wave)
% ---- amplifier ----
P.volume = 50;             % % (linear amplitude taper - assumption)
P.bass_dB = 0; P.treble_dB = 0;
P.fBass = 100;             % Hz, low-shelf corner (assumption)
P.fTreble = 10000;         % Hz, high-shelf corner (assumption)
% ---- speakers ----
P.A_out = 0.6e-3;          % m, woofer static displacement at full drive [estimate]
P.A_in = 0.15e-3;          % m, tweeter static displacement at full drive [estimate]
P.zeta_spk = 0.4;          % speaker damping ratio (total Q ~ 1.25) [estimate]
P.phi_rel_deg = 0;         % extra electrical/polarity phase of the tweeter
P.fs_out = 55;             % Hz, woofer resonance [estimate, 6.5" driver]
P.fs_in = 1800;            % Hz, tweeter resonance [estimate]
P.kappa = 0.6;             % structural coupling: rocking of tape mount by tweeter
P.r_tape = 0.05;           % m, radial distance tweeter axis -> tape
P.r_tweeter = 0.02;        % m, tweeter radiating radius
P.d_tweeter = 0.035;       % m, tweeter centre -> mirror distance
P.rho_air = 1.2;           % kg/m^3
P.c_air = 343;             % m/s
P.hd2 = 0.02;              % 2nd-harmonic distortion of cone motion (fraction) [estimate]
P.hd3 = 0.01;              % 3rd-harmonic distortion of cone motion (fraction) [estimate]
% ---- tape, spring, mirror ----
P.L = 0.0254;              % m, spring length (1 inch)
P.k = 400;                 % N/m, spring tip bending stiffness 3EI/L^3 [estimate]
P.beta_t = 1.15;           % torsion/bending stiffness ratio; helical-spring theory gives 1 + nu/2
P.m = 0.5e-3;              % kg, mirror + adhesive mass [estimate]
P.R_mirror = 6e-3;         % m, mirror radius
P.zeta = 0.01;             % modal damping ratio of the spring-mirror system [estimate]
P.C_tape = 3e-3;           % m/N, series compliance of tape + mirror/spring contact [estimate]
P.w_tape = 10e-3;          % m, lever arm converting contact compliance to rocking
P.e_h = 2e-3;              % m, lateral offset of mirror CoM from spring tip (horizontal tilt)
P.e_v = 2e-3;              % m, height of mirror CoM above spring axis (vertical tilt)
P.kappa_c = 0.2;           % cross-axis stiffness coupling (cone slant / helix)
P.cp_shift = 0.3;          % centre-of-pressure shift toward tweeter (x mirror radius)
P.eta_nl = 0.02;           % contact hardening: fractional rocking-stiffness rise at theta_ref [estimate]
P.theta_ref = 0.01;        % rad, reference tilt for eta_nl
P.theta0 = 5e-3;           % rad, kick size (angular velocity impulse = theta0 * omega_2)
P.tiltMax = 0.3;           % rad, optics clamp (outside the small-tilt contact model)
% ---- optics ----
P.D = 2.5;                 % m, mirror -> screen distance
P.brightness = 0.9;
P.spot = 3;                % px line width
P.viewHalf_cm = 10;        % cm, manual view half-width (projection scale)
P.autoFit = true;
P.exposure = 0.040;        % s, simulated exposure / visual persistence
P.noise = 0.01;            % m/s^2 RMS random base vibration (room, table)
P.noiseCorner = 300;       % Hz, corner of the first-order (Ornstein-Uhlenbeck) noise spectrum
% ---- playback ----
P.speed = 1;               % simulated seconds per wall-clock second (1 = real time, like the eye sees it)
P.sweepRate = 1;           % semitones per simulated second
P.mode = 'Physics';
P.show3D = true;
% ---- idealized Lissajous ----
P.Lis = struct('a', 3, 'b', 2, 'delta_deg', 90, 'A_cm', 8, 'B_cm', 8, ...
    'detune', 0, 'bassTreble', false);
% ---- geometry (m): x right, y out of the speaker (toward the laser), z up ----
P.geo.speakerCenter = [0; 0; 0.165];
P.geo.mirrorCenter  = [-0.035; 0.0; 0.165];   % left of the tweeter
P.geo.laserPos      = [0.12; 0.26; 0.045];    % front-right corner, elevated
P.geo.n0            = [0; 1; 0];              % mirror normal at rest
end

function St = buildStructure(P)
% Spring-mirror structural model, q = [z; theta_v; theta_h].
%  z       : spring-tip translation along the speaker axis (relative to cone)
%  theta_v : rotation about the horizontal x axis (spring torsion) -> vertical spot motion
%  theta_h : rotation about the vertical z axis (spring bending)  -> horizontal spot motion
m = P.m; L = P.L; k = P.k;
Ic = m * P.R_mirror^2 / 4;                 % disc inertia about an in-plane diameter
b = [1; -P.e_v; P.e_h];                    % y-velocity of the CoM per unit q-rate
Mm = m * (b * b.') + diag([0, Ic, Ic]);    % exact rigid-body mass matrix (small angles)

kTh = k * L^2 / 3;                         % EI/L  (k = 3EI/L^3)
kTv = P.beta_t * kTh;                      % torsional stiffness
Ks = [4 * k,       0,     -2 * k * L; ...  % Euler-Bernoulli cantilever tip matrix (EI/L^3)[12 -6L; -6L 4L^2]
      0,           kTv,    0; ...
     -2 * k * L,   0,      4 * k * L^2 / 3];
Ks(2, 3) = P.kappa_c * sqrt(kTv * kTh);
Ks(3, 2) = Ks(2, 3);
Ct = diag([P.C_tape, P.C_tape / P.w_tape^2, P.C_tape / P.w_tape^2]);
K = (Ks \ eye(3) + Ct) \ eye(3);           % series combination spring + tape/contact
K = (K + K.') / 2;

[Phi, W2] = eig(K, Mm, 'chol');
[w2, ix] = sort(real(diag(W2)));
Phi = real(Phi(:, ix));
Phi = Phi ./ sqrt(diag(Phi.' * Mm * Phi)).';     % mass-normalise
wn = sqrt(max(w2, 0));

St.M = Mm; St.K = K; St.Ks = Ks;
St.Minv = Mm \ eye(3);
St.Phi = Phi; St.wn = wn; St.fn = wn / (2 * pi);
St.zeta = P.zeta;
St.C = Mm * Phi * diag(2 * P.zeta * wn) * Phi.' * Mm;   % modal (proportional) damping
St.r_y = [1; 0; 0];                        % influence of cone translation
St.r_psi = [L; 0; 1];                      % influence of local cone rocking about z
St.r_p = [1; -P.e_v; P.e_h + P.cp_shift * P.R_mirror];  % acoustic force point
St.keff1 = 1 / (1 / k + P.C_tape);         % 1-DOF reference stiffness
St.fn1 = sqrt(St.keff1 / m) / (2 * pi);    % f_n = (1/2pi) sqrt(k/m)
% isotropic contact hardening: tau = -k3 |theta|^2 theta, k3 = eta*k_rock/theta_ref^2
St.kRock = 0.5 * (K(2, 2) + K(3, 3));
St.k3 = P.eta_nl * St.kRock / P.theta_ref^2;
end

function H = speakerFRF(f, fs, zeta)
r = f / fs;
H = 1 ./ (1 - r.^2 + 2i * zeta * r);
end

function H = eqResponse(f, P)
% First-order shelving equaliser (bass low-shelf x treble high-shelf).
s = 2i * pi * f;
GL = 10^(P.bass_dB / 20); GH = 10^(P.treble_dB / 20);
wL = 2 * pi * P.fBass; wH = 2 * pi * P.fTreble;
H = (s + wL * sqrt(GL)) ./ (s + wL / sqrt(GL)) .* ...
    GH .* (s + wH / sqrt(GH)) ./ (s + wH * sqrt(GH));
end

function F = forceFromCone(St, P, fk, So, Si)
% Generalised force phasors on q from cone-displacement phasors So, Si.
w = 2 * pi * fk;
aOut = -w.^2 .* So;                                       % woofer acceleration at the tape
aIn  = -w.^2 .* Si;                                       % tweeter acceleration
psiAcc = P.kappa * aIn / P.r_tape;                        % local rocking acceleration (rad/s^2)
pAc = P.rho_air * pi * P.r_tweeter^2 * aIn / (2 * pi * P.d_tweeter) ...
    .* exp(-1i * w * P.d_tweeter / P.c_air);              % near-field pressure with travel delay (Pa)
Fac = pi * P.R_mirror^2 * pAc;                            % force on mirror face (N)
F = -(St.M * St.r_y) * aOut - (St.M * St.r_psi) * psiAcc - St.r_p * Fac;
end

function Q = modalSolve(St, F, fk)
w = 2 * pi * fk;
den = St.wn.^2 - w.^2 + 2i * St.zeta * St.wn .* w;        % modal receptance denominators
Q = St.Phi * ((St.Phi.' * F) ./ den);
end

function [Q, So, Si, F] = harmonicResponse(St, P, f, U)
% Linear steady-state complex amplitudes for drive phasors U at frequencies f.
So = P.A_out * speakerFRF(f, P.fs_out, P.zeta_spk) .* U;
Si = P.A_in * speakerFRF(f, P.fs_in, P.zeta_spk) .* U * exp(1i * P.phi_rel_deg * pi / 180);
F = forceFromCone(St, P, f, So, Si);
Q = modalSolve(St, F, f);
end

function [n, b] = waveformHarmonics(wave, f, fBW, nMax)
% Fourier sine-series coefficients: u(t) = sum b_n sin(2 pi n f t).
N = max(1, min(nMax, floor(fBW / f)));
n = 1:N;
odd = mod(n, 2) == 1;
switch lower(wave)
    case 'sine'
        b = double(n == 1);
    case 'square'
        b = (4 / pi) ./ n .* odd;
    case 'triangle'
        sgn = ones(size(n)); sgn(mod(n, 4) == 3) = -1;
        b = (8 / pi^2) ./ n.^2 .* odd .* sgn;
    case 'sawtooth'
        b = (2 / pi) * (-1).^(n + 1) ./ n;
    otherwise
        error('Unknown waveform %s', wave);
end
keep = b ~= 0;
n = n(keep); b = b(keep);
end

function [n, So, Si] = addDistortion(n, So, Si, h2, h3)
% Weak cone nonlinearity: a quadratic / cubic term acting on the
% fundamental cone motion S1 adds h2*S1^2/|S1| at 2f and h3*S1^3/|S1|^2
% at 3f (magnitudes h*|S1|, phases 2 and 3 times the fundamental phase).
i1 = find(n == 1, 1);
if isempty(i1), return; end
orders = [2 3]; hs = [h2 h3];
S1o = So(i1); S1i = Si(i1);
for j = 1:2
    if hs(j) == 0, continue; end
    o = orders(j);
    idx = find(n == o, 1);
    if isempty(idx)
        n(end + 1) = o; So(end + 1) = 0; Si(end + 1) = 0; %#ok<AGROW>
        idx = numel(n);
    end
    So(idx) = So(idx) + hs(j) * S1o^o / max(abs(S1o), realmin)^(o - 1);
    Si(idx) = Si(idx) + hs(j) * S1i^o / max(abs(S1i), realmin)^(o - 1);
end
end

function Mdl = buildModel(P, St)
% Excitation tables for the current tone and the integrator matrices.
[n, b] = waveformHarmonics(P.wave, P.f, P.audioBandwidth, P.nHarmMax);
fk = n * P.f;
U = (P.volume / 100) * eqResponse(fk, P) .* (-1i * b);     % Re(-i b e^{iwt}) = b sin(wt)
So = P.A_out * speakerFRF(fk, P.fs_out, P.zeta_spk) .* U;
Si = P.A_in * speakerFRF(fk, P.fs_in, P.zeta_spk) .* U * exp(1i * P.phi_rel_deg * pi / 180);
[n, So, Si] = addDistortion(n, So, Si, P.hd2, P.hd3);
fk = n * P.f;
F = forceFromCone(St, P, fk, So, Si);
Q = modalSolve(St, F, fk);
% keep harmonics that matter mechanically (drive table) or visibly (speaker plot)
qm = max(abs(Q), [], 1); sm = max(abs([So; Si]), [], 1);
if max(qm) > 0, keepF = qm >= 1e-4 * max(qm); else, keepF = false(size(n)); end
if max(sm) > 0, keepS = sm >= 1e-3 * max(sm); else, keepS = false(size(n)); end
keepF(n == 1) = true; keepS(n == 1) = true;
nmax = max(n(keepF | keepS));
Np = min(8192, max(512, 32 * nmax));
keepF = keepF & n <= Np / 16; keepS = keepS & n <= Np / 16;   % table resolves >= 16 pts/period
Mdl.f = P.f;
Mdl.nF = n(keepF); Mdl.fkF = fk(keepF); Mdl.QF = Q(:, keepF); Mdl.FF = F(:, keepF);
Mdl.nS = n(keepS); Mdl.SoS = So(keepS); Mdl.SiS = Si(keepS);
Mdl.So1 = So(n == 1); Mdl.Si1 = Si(n == 1);
Mdl.nmaxF = max(Mdl.nF);
Mdl.Np = Np;
Mdl.Ft = St.Minv * periodicTable(Mdl.FF, Mdl.nF, Np);   % 3 x Np acceleration-forcing table (one period)
Mdl.qL = periodicTable(Mdl.QF, Mdl.nF, Np);              % 3 x Np small-signal steady state (linear)
Mdl.sL = periodicTable([Mdl.SoS; Mdl.SiS], Mdl.nS, Np);  % 2 x Np speaker displacements
Mdl.A = St.Minv * St.K;
Mdl.B = St.Minv * St.C;
Mdl.k3 = St.k3;
Mdl.mi2 = St.Minv(:, 2); Mdl.mi3 = St.Minv(:, 3);
Mdl.fnMax = max(St.fn);
Mdl.M = St.M; Mdl.K = St.K; Mdl.C = St.C; Mdl.Minv = St.Minv;
Mdl.Phi = St.Phi; Mdl.wn = St.wn; Mdl.zeta = St.zeta;
Mdl.tr = struct('active', false, 't0', 0, 'eta0', zeros(3, 1), 'etad0', zeros(3, 1));
end

function tab = periodicTable(C, n, Np)
% One period of real(sum_j C(:,j) exp(2i*pi*n(j)*k/Np)), k = 0..Np-1, by inverse FFT.
spec = zeros(size(C, 1), Np);
spec(:, n + 1) = C;
tab = real(ifft(spec, [], 2)) * Np;
end

function v = tableEval(tab, ph)
% Periodic linear interpolation of a one-period table at phases ph (cycles).
Np = size(tab, 2);
p = mod(ph, 1) * Np;
i0 = floor(p); fr = p - i0;
i0 = mod(i0, Np) + 1;
i1 = mod(i0, Np) + 1;
v = tab(:, i0) .* (1 - fr) + tab(:, i1) .* fr;
end

function dt = stepSize(Mdl, P, amp)
% RK4 step resolving the drive (48 pts/period), its highest retained
% harmonic (16 pts/period) and the fastest mode (30 pts/period). The
% hardening contact stiffens with amplitude, so the mode limit rises too.
fFast = Mdl.fnMax * sqrt(1 + 3 * P.eta_nl * (amp / P.theta_ref)^2);
dt = min([1 / (48 * Mdl.f), 1 / (16 * Mdl.f * Mdl.nmaxF), 1 / (30 * fFast)]);
end

function [x, ph, an, rec] = integrateChunk(x, ph, an, nSteps, dt, Mdl, rho, xi)
% Fixed-step RK4 of  q'' = Minv*F(phase) - a_n*r_y - A q - B q' - Minv*f_nl(q)
% written with scalars for speed. The tone phase (cycles) advances by f*dt
% per step, so changing f between chunks never makes the drive jump.
% a_n is a first-order (Ornstein-Uhlenbeck) random base acceleration, held
% constant over each step. rec = [phase; z; theta_v; theta_h] after each step.
Ft = Mdl.Ft; Np = size(Ft, 2);
A = Mdl.A; B = Mdl.B;
a11 = A(1, 1); a12 = A(1, 2); a13 = A(1, 3);
a21 = A(2, 1); a22 = A(2, 2); a23 = A(2, 3);
a31 = A(3, 1); a32 = A(3, 2); a33 = A(3, 3);
b11 = B(1, 1); b12 = B(1, 2); b13 = B(1, 3);
b21 = B(2, 1); b22 = B(2, 2); b23 = B(2, 3);
b31 = B(3, 1); b32 = B(3, 2); b33 = B(3, 3);
m12 = Mdl.mi2(1); m22 = Mdl.mi2(2); m32 = Mdl.mi2(3);
m13 = Mdl.mi3(1); m23 = Mdl.mi3(2); m33 = Mdl.mi3(3);
k3 = Mdl.k3;
dph = Mdl.f * dt;
q1 = x(1); q2 = x(2); q3 = x(3); v1 = x(4); v2 = x(5); v3 = x(6);
cst = [0 0.5 0.5 1];
wst = [1 2 2 1];
rec = zeros(4, nSteps);
for kStep = 1:nSteps
    an = rho * an + xi(kStep);
    pq1 = 0; pq2 = 0; pq3 = 0; pv1 = 0; pv2 = 0; pv3 = 0;
    sq1 = 0; sq2 = 0; sq3 = 0; sv1 = 0; sv2 = 0; sv3 = 0;
    for s = 1:4
        h = cst(s) * dt;
        y1 = q1 + h * pq1; y2 = q2 + h * pq2; y3 = q3 + h * pq3;
        u1 = v1 + h * pv1; u2 = v2 + h * pv2; u3 = v3 + h * pv3;
        p = ph + cst(s) * dph;
        p = (p - floor(p)) * Np;
        i0 = floor(p); fr = p - i0; i0 = i0 + 1;
        if i0 > Np, i0 = 1; end
        i1 = i0 + 1;
        if i1 > Np, i1 = 1; end
        f1 = Ft(1, i0) + fr * (Ft(1, i1) - Ft(1, i0));
        f2 = Ft(2, i0) + fr * (Ft(2, i1) - Ft(2, i0));
        f3 = Ft(3, i0) + fr * (Ft(3, i1) - Ft(3, i0));
        g = k3 * (y2 * y2 + y3 * y3);
        n2 = g * y2; n3 = g * y3;
        c1 = f1 - an - (a11 * y1 + a12 * y2 + a13 * y3) - (b11 * u1 + b12 * u2 + b13 * u3) - (m12 * n2 + m13 * n3);
        c2 = f2 - (a21 * y1 + a22 * y2 + a23 * y3) - (b21 * u1 + b22 * u2 + b23 * u3) - (m22 * n2 + m23 * n3);
        c3 = f3 - (a31 * y1 + a32 * y2 + a33 * y3) - (b31 * u1 + b32 * u2 + b33 * u3) - (m32 * n2 + m33 * n3);
        pq1 = u1; pq2 = u2; pq3 = u3;
        pv1 = c1; pv2 = c2; pv3 = c3;
        w = wst(s);
        sq1 = sq1 + w * u1; sq2 = sq2 + w * u2; sq3 = sq3 + w * u3;
        sv1 = sv1 + w * c1; sv2 = sv2 + w * c2; sv3 = sv3 + w * c3;
    end
    q1 = q1 + dt / 6 * sq1; q2 = q2 + dt / 6 * sq2; q3 = q3 + dt / 6 * sq3;
    v1 = v1 + dt / 6 * sv1; v2 = v2 + dt / 6 * sv2; v3 = v3 + dt / 6 * sv3;
    ph = ph + dph;
    ph = ph - floor(ph);
    rec(1, kStep) = ph; rec(2, kStep) = q1; rec(3, kStep) = q2; rec(4, kStep) = q3;
end
x = [q1; q2; q3; v1; v2; v3];
end

function x = ssStatePhase(Mdl, ph)
% Exact small-signal (linear) steady state [q; qdot] at tone phase ph.
E = exp(2i * pi * Mdl.nF(:) * ph);
q = real(Mdl.QF * E);
qd = real((Mdl.QF .* (2i * pi * Mdl.fkF)) * E);
x = [q; qd];
end

function x = ssStateVec(Mdl, t)
x = ssStatePhase(Mdl, mod(Mdl.f * t, 1));
end

function Mdl = setTransient(Mdl, t0, x0)
% Linear free response that makes the total state equal x0 = [q; qdot] at t0.
Mdl.tr.t0 = t0;
Mdl.tr.eta0 = Mdl.Phi.' * Mdl.M * x0(1:3);
Mdl.tr.etad0 = Mdl.Phi.' * Mdl.M * x0(4:6);
Mdl.tr.active = any(abs(x0) > 0);
end

function [q, qd] = transientModal(Mdl, tau)
% Exact damped free vibration of each mode (underdamped, zeta < 1).
w = Mdl.wn; z = Mdl.zeta;
wd = w * sqrt(1 - z^2);
tau = max(tau, 0);
Ex = exp(-z * w .* tau); Cs = cos(wd .* tau); Sn = sin(wd .* tau);
e0 = Mdl.tr.eta0; ed0 = Mdl.tr.etad0;
eta = Ex .* (e0 .* Cs + (ed0 + z * w .* e0) ./ wd .* Sn);
q = Mdl.Phi * eta;
if nargout > 1
    etad = Ex .* (ed0 .* Cs - (w.^2 .* e0 + z * w .* ed0) ./ wd .* Sn);
    qd = Mdl.Phi * etad;
end
end

function x = segmentState(Mdl, t)
% Exact linear solution: steady state + modal transient.
x = ssStateVec(Mdl, t);
if Mdl.tr.active
    [q, qd] = transientModal(Mdl, t - Mdl.tr.t0);
    x = x + [q; qd];
end
end

function [D, truncated] = selectSamples(D, NsMax, fFast)
% Reduce a chronological sample block to at most NsMax columns without
% aliasing: decimate only while >= 12 samples per period of the fastest
% content remain; otherwise keep just the most recent NsMax samples.
truncated = false;
n = size(D, 2);
if n <= NsMax, return; end
dtr = (D(1, end) - D(1, 1)) / (n - 1);
spp = 1 / max(dtr * fFast, eps);
d = min(max(1, floor(spp / 12)), ceil(n / NsMax));
D = D(:, n:-d:1);
D = D(:, end:-1:1);
if size(D, 2) > NsMax
    D = D(:, end - NsMax + 1:end);
    truncated = true;
end
end

function tv = sampleTimes(tNow, W, dt, NsMax)
% Rolling window ending at tNow (used by the idealized Lissajous view).
Ns = floor(W / dt) + 1;
Ns = max(2, min(Ns, NsMax));
tv = tNow - (Ns - 1:-1:0) * dt;
end

function [X, Y, tiltWarn] = qToScreen(q, G, P)
thV = q(2, :); thH = q(3, :);
tiltWarn = any(abs([thV thH]) > 0.25);
thV = min(max(thV, -P.tiltMax), P.tiltMax);        % optics clamp only
thH = min(max(thH, -P.tiltMax), P.tiltMax);
[X, Y] = laserSpot(G, q(1, :), thV, thH);
X = 100 * X; Y = 100 * Y;                          % m -> cm
end

function [X, Y, tv] = lissajousTrail(tNow, P, NsMax)
% Idealized x = A sin(a w t + delta), y = B sin(b w t)  (NOT the physics).
L = P.Lis;
A = L.A_cm; B = L.B_cm;
if L.bassTreble       % hypothetical: x driven by "bass", y by "treble"
    A = A * 10^(P.bass_dB / 20);
    B = B * 10^(P.treble_dB / 20);
end
fx = L.a * P.f + L.detune; fy = L.b * P.f;
dt = 1 / (96 * max(fx, fy));
tv = sampleTimes(tNow, P.exposure, dt, NsMax);
X = A * sin(2 * pi * fx * tv + L.delta_deg * pi / 180);
Y = B * sin(2 * pi * fy * tv);
end

function G = buildGeometry(P)
G.M0 = P.geo.mirrorCenter;
G.L0 = P.geo.laserPos;
G.n0 = P.geo.n0 / norm(P.geo.n0);
G.din = (G.M0 - G.L0) / norm(G.M0 - G.L0);
G.dout0 = reflectDir(G.din, G.n0);
G.w = G.dout0;                                  % screen normal (faces the mirror)
G.u = cross(G.w, [0; 0; 1]); G.u = G.u / norm(G.u);   % screen horizontal (viewer's right)
G.v = cross(G.u, G.w);                          % screen vertical (up)
G.D = P.D;
G.S0 = G.M0 + P.D * G.w;                        % nominal spot (screen origin)
G.e_h = P.e_h; G.e_v = P.e_v;
end

function dout = reflectDir(din, n)
% Vector law of reflection, column vectors (implicit expansion).
dout = din - 2 * sum(din .* n, 1) .* n;
end

function n = rotateNormal(n0, thV, thH)
% n = Rz(theta_h) * Rx(theta_v) * n0, vectorised over samples.
cv = cos(thV); sv = sin(thV); ch = cos(thH); sh = sin(thH);
x1 = n0(1) * ones(size(thV));
y1 = n0(2) * cv - n0(3) * sv;
z1 = n0(2) * sv + n0(3) * cv;
n = [x1 .* ch - y1 .* sh; x1 .* sh + y1 .* ch; z1];
end

function [X, Y, valid, Hit, dout, n, Pscr] = laserSpot(G, z, thV, thH)
% Exact ray trace: laser -> moving/tilting mirror plane -> fixed screen plane.
n = rotateNormal(G.n0, thV, thH);
dy = z + G.e_h .* thH - G.e_v .* thV;           % mirror-centre displacement along y
Pc = G.M0 + [zeros(size(dy)); dy; zeros(size(dy))];
denom = sum(G.din .* n, 1);                      % d_in . n
s = sum((Pc - G.L0) .* n, 1) ./ denom;
Hit = G.L0 + G.din .* s;
dout = G.din - 2 * denom .* n;
cosw = sum(dout .* G.w, 1);
valid = denom < -1e-3 & cosw > 0.05 & isfinite(s);
tt = sum((G.S0 - Hit) .* G.w, 1) ./ cosw;
Pscr = Hit + dout .* tt;
rel = Pscr - G.S0;
X = sum(rel .* G.u, 1);
Y = sum(rel .* G.v, 1);
X(~valid) = NaN; Y(~valid) = NaN;
end

function [e1, e2] = planeBasis(n)
n = n / norm(n);
if abs(n(3)) < 0.9, a = [0; 0; 1]; else, a = [1; 0; 0]; end
e1 = cross(n, a); e1 = e1 / norm(e1);
e2 = cross(n, e1);
end

function y = coneY(r)
% Cone surface depth (m) vs radius: 25 mm deep at the tweeter, 0 at the rim.
y = -0.025 * (0.075 - r) / (0.075 - 0.022);
end

function st = patternStats(X, Y)
ok = isfinite(X) & isfinite(Y);
st = struct('len', 0, 'wid', 0, 'ang', 0);
if nnz(ok) < 3, return; end
x = X(ok) - mean(X(ok)); y = Y(ok) - mean(Y(ok));
Cv = [mean(x.^2) mean(x.*y); mean(x.*y) mean(y.^2)];
[V, Dg] = eig(Cv);
[lam, ix] = sort(diag(Dg), 'descend');
V = V(:, ix);
lam = max(lam, 0);
st.len = 2 * sqrt(2 * lam(1));      % full extent for sinusoidal motion
st.wid = 2 * sqrt(2 * lam(2));
a = atan2(V(2, 1), V(1, 1)) * 180 / pi;
if a > 90, a = a - 180; elseif a < -90, a = a + 180; end
st.ang = a;
end

% ---------------- musical notes -----------------------------------------
function f = noteFrequency(midi)
f = 440 * 2.^((midi - 69) / 12);
end

function [name, cents, midi] = nearestNote(f)
mExact = 69 + 12 * log2(f / 440);
midi = round(mExact);
cents = 100 * (mExact - midi);
name = noteName(midi);
end

function s = noteName(midi)
names = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};
s = sprintf('%s%d', names{mod(midi, 12) + 1}, floor(midi / 12) - 1);
end

function [items, data] = noteItems()
mids = 21:119;                       % A0 ... B8
items = cell(1, numel(mids) + 1);
data = [0, mids];
items{1} = '(custom frequency)';
for i = 1:numel(mids)
    items{i + 1} = sprintf('%s  (%.2f Hz)', noteName(mids(i)), noteFrequency(mids(i)));
end
end

function s = formatFreq(f)
s = regexprep(sprintf('%.3f', f), '\.?0+$', '');
end

function s = passFail(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function C = themeColors()
C.bg = [0.08 0.085 0.10];
C.panel = [0.13 0.14 0.165];
C.field = [0.19 0.20 0.235];
C.button = [0.24 0.25 0.30];
C.text = [0.90 0.90 0.92];
C.value = [0.55 0.85 1.00];
C.accent = [1.00 0.42 0.32];
C.warn = [1.00 0.78 0.30];
C.screen = [0.018 0.018 0.022];
end

function g = greekChars()
g.theta = char(952); g.zeta = char(950); g.kappa = char(954); g.delta = char(948);
g.phi = char(966); g.mu = char(181); g.deg = char(176); g.times = char(215); g.pi = char(960);
g.beta = char(946); g.eta = char(951);
end

% =====================================================================
%                        NUMERICAL VERIFICATION
% =====================================================================
function R = runVerification(P)
R = struct('id', {}, 'name', {}, 'pass', {}, 'metric', {}, 'detail', {});
P.noise = 0;

% 1. stiffness -> natural frequency ------------------------------------
Pa = P; Pa.C_tape = 0;
Sa = buildStructure(Pa);
Pb = Pa; Pb.k = 2 * Pa.k; Sb = buildStructure(Pb);
r1 = Sb.fn1 / Sa.fn1; rM = Sb.fn ./ Sa.fn;
err = max(abs([r1; rM] - sqrt(2)));
Pc = P; Pc.k = 1.5 * P.k; Sc = buildStructure(Pc); S0 = buildStructure(P);
R(end + 1) = mkres(1, 'Higher spring stiffness raises natural frequency (m fixed)', ...
    err < 1e-7 && all(Sc.fn > S0.fn) && Sc.fn1 > S0.fn1, err, ...
    sprintf('k x2 (no tape): f_n %.2f -> %.2f Hz, ratio %.9f (expect sqrt 2); coupled ratios %s; with tape k x1.5 all modes rise', ...
    Sa.fn1, Sb.fn1, r1, mat2str(rM.', 6)));

% 2. mass -> natural frequency ----------------------------------------
Pm = P; Pm.m = 2 * P.m; Sm = buildStructure(Pm);
rm = [Sm.fn1 / S0.fn1; Sm.fn ./ S0.fn];
err = max(abs(rm - 1 / sqrt(2)));
R(end + 1) = mkres(2, 'Higher mirror mass lowers natural frequency', ...
    err < 1e-7, err, sprintf('m x2: 1-DOF and all coupled ratios %s (expect 0.7071068)', mat2str(rm.', 7)));

% 3. damping reduces the resonant peak --------------------------------
zs = [0.01 0.03 0.1];
fs = logspace(log10(10), log10(3000), 6000);
pk = zeros(size(zs));
for i = 1:numel(zs)
    Pz = P; Pz.zeta = zs(i); Sz = buildStructure(Pz);
    Q = harmonicResponse(Sz, Pz, fs, ones(size(fs)));
    pk(i) = max(abs(Q(2, :)));
end
zeta1 = 0.05; r = linspace(0.5, 1.5, 200001);
Hb = r.^2 ./ sqrt((1 - r.^2).^2 + (2 * zeta1 * r).^2);
peakTheory = 1 / (2 * zeta1 * sqrt(1 - zeta1^2));
err = abs(max(Hb) - peakTheory) / peakTheory;
R(end + 1) = mkres(3, 'More damping lowers the resonant peak', ...
    all(diff(pk) < 0) && err < 1e-6, err, ...
    sprintf('3-DOF peak |theta_y| for zeta=%s: %s mrad; 1-DOF base-excited peak %.6f vs 1/(2 zeta sqrt(1-zeta^2)) = %.6f', ...
    mat2str(zs), mat2str(1e3 * pk, 4), max(Hb), peakTheory));

% 4. phase lag vs driving frequency -----------------------------------
rr = [0.1 1 10];
ph1 = angle(1 ./ (1 - rr.^2 + 2i * zeta1 * rr)) * 180 / pi;
Fz = -(S0.M * S0.r_y);
w = 2 * pi * [0.3 1 3] * S0.fn(1);
Qz = S0.Phi * ((S0.Phi.' * Fz) ./ (S0.wn.^2 - w.^2 + 2i * S0.zeta * S0.wn .* w));
phz = unwrap(angle(Qz(1, :) ./ Fz(1))) * 180 / pi;
err = abs(ph1(2) + 90);
R(end + 1) = mkres(4, 'Mirror phase lag changes with driving frequency', ...
    err < 1e-9 && all(diff(ph1) < 0) && abs(phz(3) - phz(1)) > 90, err, ...
    sprintf('1-DOF phase at r=0.1,1,10: %s deg (exactly -90 at r=1); 3-DOF z phase at 0.3,1,3 x f1: %s deg', ...
    mat2str(ph1, 4), mat2str(phz, 4)));

% 5. law of reflection ----------------------------------------------
rs = RandStream('mt19937ar', 'Seed', 1);
N = 2000;
din = randn(rs, 3, N); din = din ./ vecnorm(din);
nn = randn(rs, 3, N); nn = nn ./ vecnorm(nn);
flip = sum(din .* nn, 1) > 0; nn(:, flip) = -nn(:, flip);
dout = reflectDir(din, nn);
eNorm = max(abs(vecnorm(dout) - 1));
eAng = max(abs(sum(dout .* nn, 1) + sum(din .* nn, 1)));
ePlane = max(abs(sum(cross(din, nn) .* dout, 1)));
err = max([eNorm eAng ePlane]);
R(end + 1) = mkres(5, 'Reflected direction obeys the vector reflection law', err < 1e-12, err, ...
    sprintf('%d random rays: max ||d_out|-1| = %.2g, max |d_out.n + d_in.n| = %.2g, coplanarity %.2g', ...
    N, eNorm, eAng, ePlane));

% 6. small rotation -> twice the deflection ---------------------------
G = buildGeometry(P);
ax = cross(G.din, G.n0); ax = ax / norm(ax);
dl = 1e-4;
n2 = rodrigues(G.n0, ax, dl);
d1 = reflectDir(G.din, G.n0); d2 = reflectDir(G.din, n2);
angOut = atan2(norm(cross(d1, d2)), dot(d1, d2));
ratio = angOut / dl;
R(end + 1) = mkres(6, 'Small mirror rotation deflects the beam by twice the angle', ...
    abs(ratio - 2) < 1e-6, abs(ratio - 2), ...
    sprintf('mirror rotated %.1e rad about the axis normal to the plane of incidence: beam turned %.6e rad (ratio %.8f)', ...
    dl, angOut, ratio));

% 7. screen distance -------------------------------------------------
Ds = [1 2 4]; disp7 = zeros(size(Ds));
for i = 1:numel(Ds)
    Pd = P; Pd.D = Ds(i); Gd = buildGeometry(Pd);
    [Xd, Yd] = laserSpot(Gd, 0, 1e-3, 0);
    disp7(i) = hypot(Xd, Yd);
end
err = max(abs(disp7(2:end) ./ disp7(1:end - 1) - 2));
R(end + 1) = mkres(7, 'Larger screen distance gives larger spot displacement', ...
    all(diff(disp7) > 0) && err < 0.01, err, ...
    sprintf('theta_y = 1 mrad: displacement at D = 1, 2, 4 m = %s mm (ratios ~2)', mat2str(1e3 * disp7, 5)));

% 8. equal-frequency Lissajous shapes ---------------------------------
t = (0:3999) / 4000;                % exactly one period, no duplicated end point
lis = @(A, B, a, b, d) deal(A * sin(2 * pi * a * t + d), B * sin(2 * pi * b * t));
[x, y] = lis(1, 1, 1, 1, 0);        eLine = max(abs(x - y));
[x, y] = lis(1, 1, 1, 1, pi / 2);   eCirc = max(abs(x.^2 + y.^2 - 1));
[x, y] = lis(1, 1, 1, 1, pi / 4);   ev = sort(eig(cov(x, y)));
eEll = abs(sqrt(ev(1) / ev(2)) - tan(pi / 8));
[x, y] = lis(2, 1, 1, 1, pi / 2);   eEll2 = max(abs((x / 2).^2 + y.^2 - 1));
err = max([eLine eCirc eEll eEll2]);
R(end + 1) = mkres(8, 'Equal-frequency Lissajous: line / circle / ellipse', err < 1e-6, err, ...
    sprintf('delta=0 line err %.1e; delta=90 circle err %.1e; delta=45 axis ratio err %.1e (expect tan 22.5); A=2B ellipse err %.1e', ...
    eLine, eCirc, eEll, eEll2));

% 9. rational ratios close ------------------------------------------
ratios = [3 2; 5 4; 1 3; 2 1];
eClose = 0;
for i = 1:size(ratios, 1)
    a = ratios(i, 1); b = ratios(i, 2);
    T = 1 / gcd(a, b);
    x0 = sin(0.7); y0 = 0;
    xT = sin(2 * pi * a * T + 0.7); yT = sin(2 * pi * b * T);
    eClose = max([eClose abs(xT - x0) abs(yT - y0)]);
end
tt = linspace(0, 50, 200001);
xi = sin(2 * pi * sqrt(2) * tt + 0.7); yi = sin(2 * pi * tt);
irr = min(hypot(xi(1000:end) - xi(1), yi(1000:end) - yi(1)));   % does not exactly return
R(end + 1) = mkres(9, 'Rational frequency ratios give closed Lissajous curves', eClose < 1e-9, eClose, ...
    sprintf('ratios 3:2, 5:4, 1:3, 2:1 return to the start after T = 1/gcd(a,b): max error %.1e; sqrt2:1 closest return over 50 periods %.2e (not closed)', ...
    eClose, irr));

% 10. full frequency range (nonlinear integrator, all waveforms) -------
fList = unique([1, 1.5, logspace(0, log10(P.fMax), 30), 440, P.fMax]);
waves = {'Sine', 'Square', 'Triangle', 'Sawtooth'};
nBad = 0; nTot = 0; maxTh = 0; maxNp = 0; minDt = inf;
Pn = P; Pn.noise = 0.05;
Sn = buildStructure(Pn); Gn = buildGeometry(Pn);
rsn = RandStream('mt19937ar', 'Seed', 7);
for iw = 1:numel(waves)
    for f = fList
        Pn.f = f; Pn.wave = waves{iw};
        Mt = buildModel(Pn, Sn);
        x0 = ssStatePhase(Mt, 0);
        amp = max(sum(abs(Mt.QF(2:3, :)), 2));
        dt = stepSize(Mt, Pn, amp);
        rho = exp(-2 * pi * Pn.noiseCorner * dt);
        xn = Pn.noise * sqrt(1 - rho^2) * randn(rsn, 300, 1);
        [xe, ~, ~, rec] = integrateChunk(x0, 0, 0, 300, dt, Mt, rho, xn);
        [X, Y] = qToScreen(rec(2:4, :), Gn, Pn);
        nTot = nTot + 1;
        bad = any(~isfinite(X)) || any(~isfinite(Y)) || any(~isfinite(xe)) || ...
            ~(dt > 0) || size(Mt.Ft, 2) > 8192;
        nBad = nBad + bad;
        maxTh = max(maxTh, max(abs(rec(3:4, :)), [], 'all'));
        maxNp = max(maxNp, size(Mt.Ft, 2));
        minDt = min(minDt, dt);
    end
end
fr = [1 27.5 440 20154];
eMap = max(abs(10.^(log10(fr)) - fr) ./ fr);
R(end + 1) = mkres(10, 'Frequency controls stable over 1 - 20,154 Hz (all waveforms)', ...
    nBad == 0 && eMap < 1e-12, nBad, ...
    sprintf('%d frequency/waveform cases integrated (nonlinear + noise), %d invalid; max |theta| %.3g mrad; table <= %d pts; min dt %.3g us; log-slider round-trip error %.1e', ...
    nTot, nBad, 1e3 * maxTh, maxNp, 1e6 * minDt, eMap));

% 11. time integrator vs exact linear solution and vs ode45 (nonlinear) --
Po = P; Po.f = 60; Po.wave = 'Sine'; Po.eta_nl = 0;
So = buildStructure(Po);
Mo = buildModel(Po, So);
Mo = setTransient(Mo, 0, -ssStateVec(Mo, 0));          % switched on at t=0, mirror at rest
dt = stepSize(Mo, Po, 0);
nS = round(0.12 / dt);
[~, ~, ~, rec] = integrateChunk(zeros(6, 1), 0, 0, nS, dt, Mo, 0, zeros(nS, 1));
qa = zeros(3, nS);
for i = 1:nS, xa = segmentState(Mo, i * dt); qa(:, i) = xa(1:3); end
scale = max(abs(qa), [], 2);
eLin = max(max(abs(rec(2:4, :) - qa), [], 2) ./ scale);

Pq = P; Pq.f = 380; Pq.wave = 'Sine'; Pq.eta_nl = 0.05;
Sq = buildStructure(Pq);
Mq = buildModel(Pq, Sq);
x0 = ssStatePhase(Mq, 0);
dt = stepSize(Mq, Pq, max(sum(abs(Mq.QF(2:3, :)), 2)));
nS = round(0.05 / dt);
[~, ~, ~, rec] = integrateChunk(x0, 0, 0, nS, dt, Mq, 0, zeros(nS, 1));
wk = 2 * pi * Mq.fkF(:);
FF = Mq.FF; Minv = Mq.Minv; Kq = Mq.K; Cq = Mq.C; k3 = Mq.k3;
rhs = @(tq, x) [x(4:6); Minv * (real(FF * exp(1i * wk * tq)) - Kq * x(1:3) - Cq * x(4:6) ...
    - k3 * (x(2)^2 + x(3)^2) * [0; x(2); x(3)])];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-14);
tq = (0:nS) * dt;
[~, xo] = ode45(rhs, tq, x0, opts);
qo = xo(2:end, 1:3).';
scale = max(abs(qo), [], 2);
eNl = max(max(abs(rec(2:4, :) - qo), [], 2) ./ scale);
R(end + 1) = mkres(11, 'Real-time integrator matches exact linear solution and ode45 (nonlinear)', ...
    eLin < 2e-3 && eNl < 5e-3, max(eLin, eNl), ...
    sprintf(['60 Hz sine switched on from rest (linear): max rel. error vs exact steady state + modal transient %.2e; ' ...
    '380 Hz with contact hardening eta = 0.05: max rel. error vs ode45 %.2e (forcing-table interpolation dominates)'], eLin, eNl));

% 12. 12-TET note table ---------------------------------------------
ref = [21 27.5; 60 261.6255653; 69 440; 119 7902.132820];
eN = max(abs(noteFrequency(ref(:, 1)) - ref(:, 2)) ./ ref(:, 2));
[nm, c] = nearestNote(440.5);
R(end + 1) = mkres(12, 'Musical note frequencies (A0..B8, 12-TET)', eN < 1e-9 && strcmp(nm, 'A4'), eN, ...
    sprintf('A0, C4, A4, B8 relative error %.1e; 440.5 Hz -> %s %+.2f cents', eN, nm, c));

% 13. nonlinear precession of elliptical orbits (the "3D rotation") -----
f0 = 300; w0 = 2 * pi * f0; a = 0.02; b = 0.01;
etaT = 0.01; beta = etaT * w0^2 / P.theta_ref^2;
Om = beta * a * b / (4 * w0);                         % first-order averaging theory
[omE, omZero, omLine] = deal(precessionRate(w0, beta, a, b), precessionRate(w0, 0, a, b), ...
    precessionRate(w0, beta, a, 0));
ratioP = abs(omE) / Om;
R(end + 1) = mkres(13, 'Contact hardening makes elliptical orbits precess at the theoretical rate', ...
    abs(ratioP - 1) < 0.1 && abs(omZero) < 0.01 * Om && abs(omLine) < 0.01 * Om, abs(ratioP - 1), ...
    sprintf(['isotropic tilt pair, 300 Hz, a = 20 mrad, b = 10 mrad, eta = %.2f: measured %.3f rad/s vs ' ...
    'beta*a*b/(4w) = %.3f rad/s (ratio %.3f); linear (eta = 0): %.2g rad/s; straight line (b = 0): %.2g rad/s'], ...
    etaT, abs(omE), Om, ratioP, omZero, omLine));

% 14. speaker harmonic distortion -------------------------------------
Pd = P; Pd.wave = 'Sine'; Pd.f = 200; Pd.hd2 = 0.03; Pd.hd3 = 0.01;
Md = buildModel(Pd, S0);
r2 = abs(Md.SoS(Md.nS == 2)) / abs(Md.So1);
r3 = abs(Md.SoS(Md.nS == 3)) / abs(Md.So1);
Pd.hd2 = 0; Pd.hd3 = 0;
Md0 = buildModel(Pd, S0);
if isempty(r2), r2 = NaN; end
if isempty(r3), r3 = NaN; end
err = max(abs([r2 - 0.03, r3 - 0.01]));
R(end + 1) = mkres(14, 'Harmonic distortion adds exactly the set 2f / 3f cone motion', ...
    err < 1e-12 && isequal(Md0.nS, 1) && isequal(Md0.nF, 1), err, ...
    sprintf('200 Hz sine: |S(2f)|/|S(f)| = %.6f (set 0.03), |S(3f)|/|S(f)| = %.6f (set 0.01); with 0%% distortion only the fundamental is present', ...
    r2, r3));

% 15. phase-continuous frequency change -------------------------------
P1 = P; P1.f = 440; P1.wave = 'Sine'; M1 = buildModel(P1, S0);
P2 = P; P2.f = 470; P2.wave = 'Sine'; M2 = buildModel(P2, S0);
dt = min(stepSize(M1, P1, 0), stepSize(M2, P2, 0));
[x1, p1, ~, rec1] = integrateChunk(ssStatePhase(M1, 0), 0, 0, 200, dt, M1, 0, zeros(200, 1));
[~, ~, ~, rec2] = integrateChunk(x1, p1, 0, 200, dt, M2, 0, zeros(200, 1));
u = sin(2 * pi * [rec1(1, :) rec2(1, :)]);
du = max(abs(diff(u)));
bound = 2 * pi * 470 * dt;
R(end + 1) = mkres(15, 'Tone stays phase-continuous when the frequency changes', ...
    du <= bound * 1.0001, du / bound, ...
    sprintf('440 -> 470 Hz mid-run: largest sample-to-sample tone change %.4g vs smooth-signal bound 2*pi*f*dt = %.4g', du, bound));
end

function om = precessionRate(w0, beta, a, b)
% Free isotropic tilt pair theta'' + w0^2 theta + beta |theta|^2 theta = 0,
% integrated with the production RK4 routine; returns the rotation rate of
% the orbit's apsides (rad/s).
Mf.Ft = zeros(3, 512); Mf.A = w0^2 * eye(3); Mf.B = zeros(3); Mf.k3 = beta;
Mf.mi2 = [0; 1; 0]; Mf.mi3 = [0; 0; 1]; Mf.f = w0 / (2 * pi);
dt = 2 * pi / w0 / 200;
nS = 200 * 200;
x0 = [0; a; 0; 0; 0; b * w0];
[~, ~, ~, rec] = integrateChunk(x0, 0, 0, nS, dt, Mf, 0, zeros(nS, 1));
xx = rec(3, :); yy = rec(4, :); tt = (1:nS) * dt;
rr = hypot(xx, yy);
pk = find(rr(2:end - 1) > rr(1:end - 2) & rr(2:end - 1) >= rr(3:end)) + 1;
ang = mod(atan2(yy(pk), xx(pk)), pi);
d = diff(ang); d = mod(d + pi / 2, pi) - pi / 2;
u = [ang(1), ang(1) + cumsum(d)];
c = polyfit(tt(pk), u, 1);
om = c(1);
end

function r = mkres(id, name, pass, metric, detail)
r = struct('id', id, 'name', name, 'pass', logical(pass), 'metric', double(metric), 'detail', detail);
end

function v = rodrigues(v, k, a)
v = v * cos(a) + cross(k, v) * sin(a) + k * dot(k, v) * (1 - cos(a));
end

function printVerification(R)
fprintf('\n=== Laser simulator numerical verification ===\n');
for i = 1:numel(R)
    fprintf('%2d  %-4s  %-70s  metric = %.3g\n      %s\n', R(i).id, passFail(R(i).pass), ...
        R(i).name, R(i).metric, R(i).detail);
end
fprintf('%d / %d passed (internal numerical consistency only).\n\n', sum([R.pass]), numel(R));
end
