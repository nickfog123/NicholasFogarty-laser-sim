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
%     Tone      : u(t) = sum_n b_n sin(2*pi*n*f*t), band-limited Fourier series
%                 of sine / square / triangle / sawtooth.
%     Amplifier : first-order low-shelf (bass, 100 Hz) and high-shelf (treble,
%                 10 kHz) equaliser, linear master volume.
%     Speakers  : each cone is a driven damped oscillator
%                 m_i s_i'' + c_i s_i' + k_i s_i = F_i(t)
%                 written as a displacement FRF  A_i / (1 - r^2 + 2 i zeta r).
%     Mount     : base translation of the woofer cone (y_b = s_out) plus a
%                 local cone rocking angle psi_b = kappa*s_in/r_tape caused by
%                 the tweeter (structural coupling) and near-field acoustic
%                 pressure of the tweeter acting on the mirror face.
%     Mirror    : 3-DOF rigid body on a cantilevered coil spring,
%                 q = [z; theta_v; theta_h]  (z = relative tip translation
%                 along the speaker axis, theta_v = vertical tilt (about the
%                 horizontal X axis), theta_h = horizontal tilt (about the
%                 vertical Z axis)):
%                 M q'' + C q' + K q = -M r_y yb'' - M r_psi psi_b'' + F_ac
%                 K = series(spring Euler-Bernoulli tip matrix + torsion,
%                 tape/contact compliance). Proportional (modal) damping.
%     Solution  : exact frequency-domain steady state per harmonic (stored as
%                 a one-period lookup table) + exact analytical modal
%                 transient + narrow-band random base vibration. No ODE
%                 integration in the animation loop.
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
G   = buildGeometry(P);         % optical geometry
S   = struct();                 % runtime state
S.t = 0;                        % simulated time (s)
S.playing = false;
S.segs = {};                    % piecewise-constant parameter segments
S.segStart = zeros(1, 0);       % start time of each segment (s)
S.noiseSeed = makeNoiseSeed(12);
S.needRebuild = false;
S.geoDirty = false;
S.resDirty = true;
S.schemDirty = true;
S.frame = 0;
S.NsMax = 3000;                 % adaptive cap on trail samples per frame
S.lastTic = tic;
S.viewHalf = P.viewHalf_cm;     % current displayed half-width (cm)
S.lastStats = struct('len', 0, 'wid', 0, 'ang', 0);
S.saturated = false;
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
        nRows = [8 4 8 11 7 6 7];
        left.RowHeight = num2cell(panelHeight(nRows));

        buildPanelA(makePanel(left, 1, 'A  Frequency & waveform (tone generator)', nRows(1)));
        buildPanelB(makePanel(left, 2, 'B  Stereo amplifier', nRows(2)));
        buildPanelC(makePanel(left, 3, 'C  Speaker dynamics (both speakers get the same signal)', nRows(3)));
        buildPanelD(makePanel(left, 4, 'D  Tape, 1-inch spring & mirror  [estimated constants]', nRows(4)));
        buildPanelE(makePanel(left, 5, 'E  Laser & screen geometry', nRows(5)));
        buildPanelF(makePanel(left, 6, 'F  Playback & simulation', nRows(6)));
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
        H.spkInfo = addLabel(g, 8, [1 3], '', CLR.value);
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
        addSlider(g, 8, ['Initial perturb. ' GL.theta '_0'], [0 50], P.theta0 * 1e3, @(v) sprintf('%.1f mrad', v), ...
            @(v) setParam('theta0', v * 1e-3, 'none'));
        H.fnInfo1 = addLabel(g, 9, [1 3], '', CLR.value);
        H.fnInfo2 = addLabel(g, 10, [1 3], '', CLR.value);
        lg = uigridlayout(g, [1 2]);
        lg.Layout.Row = 11; lg.Layout.Column = [1 3];
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
        addSlider(g, 6, 'Disturbance (base)', [0 2], P.noise, @(v) sprintf('%.2f m/s%s', v, char(178)), ...
            @(v) setParam('noise', v, 'dyn'));
        H.autoFit = uicheckbox(g, 'Text', 'Auto-fit screen view to pattern (sizes still reported in cm)', ...
            'Value', P.autoFit, 'FontColor', CLR.text, ...
            'ValueChangedFcn', @(src, ~) setParam('autoFit', src.Value, 'none'));
        H.autoFit.Layout.Row = 7; H.autoFit.Layout.Column = [1 3];
    end

    % ---------------- Panel F: playback -------------------------------
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
        addSlider(g, 3, 'Simulation speed', [-4 0], log10(P.speed), ...
            @(v) sprintf('%.3g%s', 10^v, GL.times), @(v) setParam('speed', 10^v, 'none'));
        H.modeDrop = addDropdown(g, 4, 'Display model', ...
            {'Physics model (speaker-spring-mirror)', 'Idealized Lissajous (comparison)'}, ...
            {'Physics', 'Lissajous'}, P.mode, @(v) setParam('mode', v, 'none'));
        H.show3D = uicheckbox(g, 'Text', 'Update 3D schematic', 'Value', P.show3D, ...
            'FontColor', CLR.text, 'ValueChangedFcn', @(src, ~) setParam('show3D', src.Value, 'none'));
        H.show3D.Layout.Row = 5; H.show3D.Layout.Column = [1 2];
        H.timeLbl = addLabel(g, 6, [1 3], '', CLR.value);
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
        g.RowHeight = {'1x', 60};
        g.BackgroundColor = CLR.bg; g.Padding = [4 4 4 4];
        H.axRes = darkAxes(g, 1, 'Steady-state response to a sine at the current amplifier settings', ...
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
        if S.needRebuild
            applyDynamics();
        end
        if S.geoDirty
            G = buildGeometry(P);
            S.geoDirty = false;
            S.resDirty = true;
        end
    end

    function applyDynamics()
        % Rebuild the model and start a new segment at the current time.
        % State continuity: the old segment's state at t_c becomes the new
        % initial condition; the difference from the new steady state is
        % carried by an exact free (transient) response.
        St = buildStructure(P);
        G = buildGeometry(P);
        newM = buildModel(P, St, S.noiseSeed);
        if isempty(S.segs)
            xOld = zeros(6, 1);
        else
            xOld = segmentState(S.segs{end}, S.t);
        end
        newM = setTransient(newM, S.t, xOld - ssStateVec(newM, S.t));
        pushSegment(newM);
        S.needRebuild = false;
        S.resDirty = true;
        updateReadouts();
    end

    function pushSegment(Mdl)
        if ~isempty(S.segStart) && S.segStart(end) >= S.t
            keep = S.segStart < S.t;
            S.segs = S.segs(keep);
            S.segStart = S.segStart(keep);
        end
        S.segs{end + 1} = Mdl;
        S.segStart(end + 1) = S.t;
        horizon = S.t - 1.0;        % simulated seconds of history kept
        while numel(S.segStart) >= 2 && S.segStart(2) <= horizon
            S.segs(1) = [];
            S.segStart(1) = [];
        end
        if numel(S.segs) > 80
            S.segs = S.segs(end - 79:end);
            S.segStart = S.segStart(end - 79:end);
        end
    end

    function resetSimulation()
        % Sound switched on at t = 0 with the mirror at rest, plus the
        % user-selected initial angular perturbation theta_0.
        S.t = 0;
        S.segs = {};
        S.segStart = zeros(1, 0);
        St = buildStructure(P);
        G = buildGeometry(P);
        M0 = buildModel(P, St, S.noiseSeed);
        xStart = [0; P.theta0; 0.7 * P.theta0; 0; 0; 0];
        M0 = setTransient(M0, 0, xStart - ssStateVec(M0, 0));
        S.segs = {M0};
        S.segStart = 0;
        S.needRebuild = false;
        S.geoDirty = false;
        S.resDirty = true;
        S.schemDirty = true;
    end

    function onReset()
        resetSimulation();
        updateReadouts();
        renderFrame(true);
    end

    function onKick()
        % Impulsive angular disturbance (e.g. a tap on the speaker).
        processPending();
        Mdl = S.segs{end};
        xNow = segmentState(Mdl, S.t);
        dx = [0; P.theta0; -0.7 * P.theta0; 0; 0; 0];
        Mdl = setTransient(Mdl, S.t, xNow + dx - ssStateVec(Mdl, S.t));
        pushSegment(Mdl);
        if ~S.playing, renderFrame(true); end
    end

    % -----------------------------------------------------------------
    % Frequency controls
    % -----------------------------------------------------------------
    function setFrequency(f, source)
        if ~isfinite(f), return; end
        f = round(f * 1000) / 1000;
        f = min(max(f, P.fMin), P.fMax);
        P.f = f;
        refreshFrequencyWidgets(source);
        S.needRebuild = true;
        S.resDirty = true;
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
            processPending();
            S.t = S.t + dtWall * P.speed;
            renderFrame(false);
        catch err
            stopPlayback();
            if isvalid(H.fig)
                uialert(H.fig, sprintf('Animation stopped: %s', err.message), 'Simulation error');
            end
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
    % Rendering
    % -----------------------------------------------------------------
    function renderFrame(forceAll)
        if nargin < 1, forceAll = false; end
        tStart = tic;
        S.frame = S.frame + 1;
        Mdl = S.segs{end};

        % ---- compute the visible trajectory (rolling exposure window) ----
        if strcmp(P.mode, 'Physics')
            dt = sampleStep(Mdl, Mdl.f);
            tv = sampleTimes(S.t, P.exposure, dt, S.NsMax);
            q = evalSegments(S.segs, S.segStart, tv, false);
            [X, Y, S.saturated] = qToScreen(q, G, P);
        else
            [X, Y, tv] = lissajousTrail(S.t, P, S.NsMax);
            S.saturated = false;
        end
        Ns = numel(X);

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
            if target > S.viewHalf || target < 0.45 * S.viewHalf
                if target > S.viewHalf
                    S.viewHalf = target;
                else
                    S.viewHalf = S.viewHalf + 0.25 * (target - S.viewHalf);
                end
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
        S.lastStats = st;
        [nm, cents] = nearestNote(P.f);
        if strcmp(P.mode, 'Physics')
            modeStr = 'PHYSICS MODEL';
        elseif P.Lis.bassTreble
            modeStr = 'HYPOTHETICAL two-frequency (NOT how bass/treble work)';
        else
            modeStr = 'IDEALIZED LISSAJOUS (not the physical system)';
        end
        H.info.String = sprintf(['%s Hz  (%s %+.0f c)  %s wave\n%s\n' ...
            'pattern %.2f x %.2f cm, major axis %+.0f deg   view %s%.3g cm'], ...
            formatFreq(P.f), nm, cents, P.wave, modeStr, st.len, st.wid, st.ang, ...
            char(177), S.viewHalf);
        if S.saturated
            H.warnTxt.String = 'Tilt > contact limit: mirror would rattle on the spring (tanh limit applied)';
        elseif Ns >= S.NsMax && strcmp(P.mode, 'Physics')
            H.warnTxt.String = sprintf('High frequency: last %.3g ms of the exposure drawn (pattern is periodic)', ...
                1e3 * (tv(end) - tv(1)));
        else
            H.warnTxt.String = '';
        end
        H.timeLbl.Text = sprintf('t_sim = %.6f s   frame %d   samples %d', S.t, S.frame, Ns);
        H.status.Text = strsplit(statusText(Mdl, st), newline);

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

        % ---- performance safeguard: adapt sample budget ------------------
        el = toc(tStart);
        if el > 0.030
            S.NsMax = max(600, round(S.NsMax * 0.85));
        elseif el < 0.012
            S.NsMax = min(3000, round(S.NsMax * 1.05) + 1);
        end
    end

    function txt = statusText(Mdl, st)
        thV = abs(Mdl.Q(2, 1)); thH = abs(Mdl.Q(3, 1));
        dphi = angle(Mdl.Q(2, 1) / (Mdl.Q(3, 1) + eps)) * 180 / pi;
        txt = sprintf([' Fundamental: |%s_x| = %.3g mrad, |%s_y| = %.3g mrad, phase(%s_y - %s_x) = %+.0f%s' ...
            '   |z| = %.3g %sm\n Trail: %.2f x %.2f cm  (axis ratio %.2f)   harmonics used: %d   transient: %s'], ...
            GL.theta, 1e3 * thH, GL.theta, 1e3 * thV, GL.theta, GL.theta, dphi, GL.deg, ...
            1e6 * abs(Mdl.Q(1, 1)), GL.mu, st.len, st.wid, st.wid / max(st.len, eps), ...
            numel(Mdl.n), onOff(transientAlive(Mdl, S.t)));
    end

    function updateSignalPlots()
        Mdl = S.segs{end};
        dt = sampleStep(Mdl, Mdl.f);
        Wp = min(3 / Mdl.f, 0.5);
        Ns = min(2000, max(50, ceil(Wp / dt) + 1));
        tv = S.t - (Ns - 1:-1:0) * dt;
        [q, sp] = evalSegments(S.segs, S.segStart, tv, true);
        tm = (tv - S.t) * 1e3;
        set(H.lSpkOut, 'XData', tm, 'YData', sp(1, :) * 1e6);
        set(H.lSpkIn,  'XData', tm, 'YData', sp(2, :) * 1e6);
        set(H.lSpr,    'XData', tm, 'YData', q(1, :) * 1e6);
        set(H.lThX,    'XData', tm, 'YData', q(3, :) * 1e3);
        set(H.lThY,    'XData', tm, 'YData', q(2, :) * 1e3);
        xl = [tm(1) 0];
        if xl(2) <= xl(1), xl(1) = xl(2) - 1e-6; end
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
        [~, ix] = max(thx + thy);
        H.resNote.Text = sprintf([' Coupled natural frequencies: %s Hz  (undamped)   |   1-DOF reference ' ...
            'f_n = sqrt(k_eff/m)/2%s = %.1f Hz   |   largest tilt response at %.1f Hz. ' ...
            'Above the speaker resonances the cone acceleration is roughly constant, so the mirror ' ...
            'response falls as 1/f^2 above its own modes.'], ...
            strjoin(arrayfun(@(x) sprintf('%.1f', x), St.fn(:)', 'UniformOutput', false), ', '), ...
            GL.pi, St.fn1, fs(ix));
    end

    function updateReadouts()
        if ~isfield(H, 'eqInfo'), return; end
        Mdl = S.segs{end};
        Heq = eqResponse(P.f, P);
        H.eqInfo.Text = sprintf('EQ gain at %s Hz: %+.2f dB   (linear volume %.2f)', ...
            formatFreq(P.f), 20 * log10(abs(Heq)), P.volume / 100);
        H.spkInfo.Text = sprintf('At f: |s_out| = %.3g %sm, |s_in| = %.3g %sm (fundamental)', ...
            1e6 * abs(Mdl.So(1)), GL.mu, 1e6 * abs(Mdl.Si(1)), GL.mu);
        H.fnInfo1.Text = sprintf('1-DOF f_n = %.1f Hz (k_eff = %.0f N/m, m = %.2f g)', ...
            St.fn1, St.keff1, P.m * 1e3);
        H.fnInfo2.Text = sprintf('Coupled modes: %.1f, %.1f, %.1f Hz', St.fn(1), St.fn(2), St.fn(3));
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
            q = evalSegments(S.segs, S.segStart, S.t, false);
        else
            q = zeros(3, 1);
        end
        thV = P.theta_lim * tanh(q(2) / P.theta_lim) * exag;
        thH = P.theta_lim * tanh(q(3) / P.theta_lim) * exag;
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
                'Position', [120 120 1180 520], 'Color', CLR.bg);
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
P.nHarmMax = 300;          % max number of Fourier harmonics considered
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
% ---- tape, spring, mirror ----
P.L = 0.0254;              % m, spring length (1 inch)
P.k = 400;                 % N/m, spring tip bending stiffness 3EI/L^3 [estimate]
P.beta_t = 0.3;            % torsional/bending rotational stiffness ratio [estimate]
P.m = 0.5e-3;              % kg, mirror + adhesive mass [estimate]
P.R_mirror = 6e-3;         % m, mirror radius
P.zeta = 0.05;             % modal damping ratio of the spring-mirror system [estimate]
P.C_tape = 0.5e-3;         % m/N, series compliance of tape + mirror/spring contact
P.w_tape = 10e-3;          % m, lever arm converting tape compliance to rotation
P.e_h = 3e-3;              % m, lateral offset of mirror CoM from spring tip (horizontal tilt)
P.e_v = 6e-3;              % m, height of mirror CoM above spring axis (vertical tilt)
P.kappa_c = 0.2;           % cross-axis stiffness coupling (cone slant / helix)
P.cp_shift = 0.3;          % centre-of-pressure shift toward tweeter (x mirror radius)
P.theta0 = 5e-3;           % rad, initial angular perturbation at reset / kick
P.theta_lim = 0.15;        % rad, contact limit (mirror would lift / rattle)
% ---- optics ----
P.D = 2.5;                 % m, mirror -> screen distance
P.brightness = 0.9;
P.spot = 3;                % px line width
P.viewHalf_cm = 10;        % cm, manual view half-width (projection scale)
P.autoFit = true;
P.exposure = 0.040;        % s, simulated exposure / visual persistence
P.noise = 0.05;            % m/s^2 RMS narrow-band base vibration (room, table)
% ---- playback ----
P.speed = 0.01;            % simulated seconds per wall-clock second (slow motion)
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
St.Phi = Phi; St.wn = wn; St.fn = wn / (2 * pi);
St.zeta = P.zeta;
St.C = Mm * Phi * diag(2 * P.zeta * wn) * Phi.' * Mm;   % modal (proportional) damping
St.r_y = [1; 0; 0];                        % influence of cone translation
St.r_psi = [L; 0; 1];                      % influence of local cone rocking about z
St.r_p = [1; -P.e_v; P.e_h + P.cp_shift * P.R_mirror];  % acoustic force point
St.keff1 = 1 / (1 / k + P.C_tape);         % 1-DOF reference stiffness
St.fn1 = sqrt(St.keff1 / m) / (2 * pi);    % f_n = (1/2pi) sqrt(k/m)
end

function [Q, So, Si, F] = harmonicResponse(St, P, f, U)
% Steady-state complex amplitudes for drive phasors U at frequencies f.
w = 2 * pi * f;
Ho = speakerFRF(f, P.fs_out, P.zeta_spk);
Hi = speakerFRF(f, P.fs_in, P.zeta_spk);
So = P.A_out * Ho .* U;                                   % woofer cone displacement (m)
Si = P.A_in * Hi .* U * exp(1i * P.phi_rel_deg * pi / 180); % tweeter cone displacement (m)
aOut = -w.^2 .* So;                                       % cone acceleration
aIn  = -w.^2 .* Si;
psiAcc = P.kappa * aIn / P.r_tape;                        % local rocking acceleration (rad/s^2)
pAc = P.rho_air * pi * P.r_tweeter^2 * aIn / (2 * pi * P.d_tweeter);   % near-field pressure (Pa)
Fac = pi * P.R_mirror^2 * pAc;                            % force on mirror face (N)
F = -(St.M * St.r_y) * aOut - (St.M * St.r_psi) * psiAcc - St.r_p * Fac;
den = St.wn.^2 - w.^2 + 2i * St.zeta * St.wn .* w;        % modal receptance denominators
Q = St.Phi * ((St.Phi.' * F) ./ den);
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

function Mdl = buildModel(P, St, seed)
% Build one parameter segment: steady state (per harmonic + one-period
% lookup table), modal data for the transient and the disturbance.
[n, b] = waveformHarmonics(P.wave, P.f, P.audioBandwidth, P.nHarmMax);
fk = n * P.f;
U = (P.volume / 100) * eqResponse(fk, P) .* (-1i * b);     % Re(-i b e^{iwt}) = b sin(wt)
[Q, So, Si, F] = harmonicResponse(St, P, fk, U);
% prune harmonics that are negligible both mechanically and acoustically
qm = max(abs(Q), [], 1); sm = max(abs([So; Si]), [], 1);
keep = (qm >= 1e-5 * max(qm)) | (sm >= 1e-5 * max(sm));
if max(qm) == 0 && max(sm) == 0, keep = false(size(n)); end
keep(1) = true;
n = n(keep); fk = fk(keep); b = b(keep); U = U(keep);
Q = Q(:, keep); So = So(keep); Si = Si(keep); F = F(:, keep);

Np = min(4096, max(256, 16 * max(n)));
ph = (0:Np - 1) / Np;
E = exp(2i * pi * n(:) * ph);              % Nh x Np
Mdl.f = P.f; Mdl.n = n; Mdl.fk = fk; Mdl.b = b; Mdl.U = U;
Mdl.Q = Q; Mdl.So = So; Mdl.Si = Si; Mdl.F = F;
Mdl.Np = Np;
Mdl.qL = real(Q * E);                      % 3 x Np one-period table
Mdl.sL = real([So; Si] * E);               % 2 x Np speaker displacements
Mdl.M = St.M; Mdl.K = St.K; Mdl.C = St.C;
Mdl.Phi = St.Phi; Mdl.wn = St.wn; Mdl.zeta = St.zeta;
Mdl.tr = struct('active', false, 't0', 0, 'eta0', zeros(3, 1), 'etad0', zeros(3, 1));
% narrow-band random base vibration: sinusoids clustered around each mode
Nn = numel(seed.u);
modeIdx = mod(0:Nn - 1, 3) + 1;
wj = St.wn(modeIdx) .* (1 + 1.5 * St.zeta * seed.u(:));
wj = max(wj, 2 * pi * 0.5);
aj = P.noise * sqrt(2 / Nn) * ones(1, Nn);
Fj = -(St.M * St.r_y) * aj;
den = St.wn.^2 - (wj.').^2 + 2i * St.zeta * St.wn .* (wj.');
Mdl.noiseQ = St.Phi * ((St.Phi.' * Fj) ./ den);
Mdl.noiseW = wj(:);
Mdl.noisePh = seed.ph(:);
Mdl.noiseOn = P.noise > 0;
end

function seed = makeNoiseSeed(N)
rs = RandStream('mt19937ar', 'Seed', 2026);
seed.u = randn(rs, N, 1);
seed.ph = 2 * pi * rand(rs, N, 1);
end

function Mdl = setTransient(Mdl, t0, x0)
% Free response that makes the total state equal x0 = [q; qdot] at t0.
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

function tf = transientAlive(Mdl, t)
tf = Mdl.tr.active && exp(-Mdl.zeta * min(Mdl.wn) * max(t - Mdl.tr.t0, 0)) > 1e-6;
end

function x = ssStateVec(Mdl, t)
% Exact steady-state state vector [q; qdot] at time t.
E = exp(2i * pi * Mdl.n(:) * mod(Mdl.f * t, 1));
q = real(Mdl.Q * E);
qd = real((Mdl.Q .* (2i * pi * Mdl.fk)) * E);
x = [q; qd];
end

function x = segmentState(Mdl, t)
x = ssStateVec(Mdl, t);
if Mdl.tr.active
    [q, qd] = transientModal(Mdl, t - Mdl.tr.t0);
    x = x + [q; qd];
end
end

function [q, sp] = evalSegments(segs, segStart, t, wantSpk)
% Mirror coordinates at sample times t using the segment active at each t.
q = zeros(3, numel(t));
sp = zeros(2, numel(t));
idx = sum(segStart(:) <= t, 1);
idx(idx < 1) = 1;
for kk = unique(idx)
    mask = idx == kk;
    Mdl = segs{kk};
    tk = t(mask);
    ph = mod(Mdl.f * tk, 1) * Mdl.Np;
    i0 = floor(ph); fr = ph - i0;
    i0 = mod(i0, Mdl.Np) + 1;
    i1 = mod(i0, Mdl.Np) + 1;
    qk = Mdl.qL(:, i0) .* (1 - fr) + Mdl.qL(:, i1) .* fr;
    if Mdl.tr.active && transientAlive(Mdl, min(tk))
        qk = qk + transientModal(Mdl, tk - Mdl.tr.t0);
    end
    if Mdl.noiseOn
        qk = qk + real(Mdl.noiseQ * exp(1i * (Mdl.noiseW * tk + Mdl.noisePh)));
    end
    q(:, mask) = qk;
    if wantSpk
        sp(:, mask) = Mdl.sL(:, i0) .* (1 - fr) + Mdl.sL(:, i1) .* fr;
    end
end
end

function dt = sampleStep(Mdl, f)
% Time step that resolves the drive (64 pts/period), its highest retained
% harmonic (12 pts/period) and the mirror modes (24 pts/period).
dt = min(1 / (64 * f), 1 / (12 * Mdl.fk(end)));
if Mdl.noiseOn || Mdl.tr.active
    dt = min(dt, 1 / (24 * max(Mdl.wn) / (2 * pi)));
end
end

function tv = sampleTimes(tNow, W, dt, NsMax)
% Rolling exposure window ending at tNow. If the window would need more
% than NsMax samples, only the most recent NsMax*dt seconds are drawn
% (the periodic pattern looks identical; anti-aliasing is preserved).
Ns = floor(W / dt) + 1;
Ns = max(2, min(Ns, NsMax));
tv = tNow - (Ns - 1:-1:0) * dt;
end

function [X, Y, saturated] = qToScreen(q, G, P)
thV = q(2, :); thH = q(3, :);
saturated = any(abs([thV thH]) > 0.6 * P.theta_lim);
thV = P.theta_lim * tanh(thV / P.theta_lim);       % simplified contact limit
thH = P.theta_lim * tanh(thH / P.theta_lim);
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

function s = onOff(tf)
if tf, s = 'decaying'; else, s = 'negligible'; end
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

% 10. full frequency range --------------------------------------------
fList = unique([1, 1.5, logspace(0, log10(P.fMax), 30), 440, P.fMax]);
waves = {'Sine', 'Square', 'Triangle', 'Sawtooth'};
nBad = 0; nTot = 0; maxTh = 0; maxNs = 0; seed = makeNoiseSeed(12);
Pn = P; Pn.noise = 0.2;
Sn = buildStructure(Pn); Gn = buildGeometry(Pn);
for iw = 1:numel(waves)
    for f = fList
        Pn.f = f; Pn.wave = waves{iw};
        M = buildModel(Pn, Sn, seed);
        M = setTransient(M, 0, -ssStateVec(M, 0));
        dt = sampleStep(M, f);
        tv = sampleTimes(0.37, Pn.exposure, dt, 3000);
        q = evalSegments({M}, 0, tv, false);
        [X, Y] = qToScreen(q, Gn, Pn);
        nTot = nTot + 1;
        bad = any(~isfinite(X)) || any(~isfinite(Y)) || any(~isfinite(q(:))) || numel(tv) > 3000 || M.Np > 4096;
        nBad = nBad + bad;
        maxTh = max(maxTh, max(abs(q(2:3, :)), [], 'all'));
        maxNs = max(maxNs, numel(tv));
    end
end
fr = [1 27.5 440 20154];
eMap = max(abs(10.^(log10(fr)) - fr) ./ fr);
R(end + 1) = mkres(10, 'Frequency controls stable over 1 - 20,154 Hz (all waveforms)', ...
    nBad == 0 && eMap < 1e-12, nBad, ...
    sprintf('%d frequency/waveform cases, %d invalid; max |theta| %.3g mrad; max samples/frame %d; log-slider round-trip error %.1e', ...
    nTot, nBad, 1e3 * maxTh, maxNs, eMap));

% 11. analytic steady state + transient vs ode45 ----------------------
Po = P; Po.f = 60; Po.wave = 'Sine'; Po.noise = 0;
So = buildStructure(Po);
Mo = buildModel(Po, So, seed);
Mo = setTransient(Mo, 0, -ssStateVec(Mo, 0));
Fp = Mo.F; wk = 2 * pi * Mo.fk(:);
rhs = @(tq, x) [x(4:6); Mo.M \ (real(Fp * exp(1i * wk * tq)) - Mo.K * x(1:3) - Mo.C * x(4:6))];
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-13);
tq = linspace(0, 0.12, 400);
[~, xo] = ode45(rhs, tq, zeros(6, 1), opts);
xa = zeros(6, numel(tq));
for i = 1:numel(tq), xa(:, i) = segmentState(Mo, tq(i)); end
scale = max(abs(xo(:, 1:3)), [], 1).';
eOde = max(max(abs(xo(:, 1:3).' - xa(1:3, :)), [], 2) ./ scale);
tl = rand(rs, 1, 500) * 0.1;
ql = evalSegments({Mo}, 0, tl, false);
qe = zeros(3, numel(tl));
for i = 1:numel(tl), xe = segmentState(Mo, tl(i)); qe(:, i) = xe(1:3); end
eLut = max(max(abs(ql - qe), [], 2) ./ scale);
R(end + 1) = mkres(11, 'Analytical solution matches direct ODE integration (ode45)', ...
    eOde < 1e-4 && eLut < 2e-2, eOde, ...
    sprintf('60 Hz sine switched on at t=0 (mirror at rest), %d harmonic: max relative error vs ode45 %.2e; lookup-table interpolation error %.2e', ...
    numel(Mo.n), eOde, eLut));

% 12. 12-TET note table ---------------------------------------------
ref = [21 27.5; 60 261.6255653; 69 440; 119 7902.132820];
eN = max(abs(noteFrequency(ref(:, 1)) - ref(:, 2)) ./ ref(:, 2));
[nm, c] = nearestNote(440.5);
R(end + 1) = mkres(12, 'Musical note frequencies (A0..B8, 12-TET)', eN < 1e-9 && strcmp(nm, 'A4'), eN, ...
    sprintf('A0, C4, A4, B8 relative error %.1e; 440.5 Hz -> %s %+.2f cents', eN, nm, c));
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
    fprintf('%2d  %-4s  %-62s  metric = %.3g\n      %s\n', R(i).id, passFail(R(i).pass), ...
        R(i).name, R(i).metric, R(i).detail);
end
fprintf('%d / %d passed (internal numerical consistency only).\n\n', sum([R.pass]), numel(R));
end
