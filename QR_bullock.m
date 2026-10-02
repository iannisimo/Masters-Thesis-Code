function [circuit, err, b_time, s_time] = QR_bullock(d, n, IMAG, seed, psi_in)
    % Usable both ways:
    %   QR_bullock            % run as a "script" (Run button / matlab -batch QR_bullock), uses defaults
    %   [c, e, b, s] = QR_bullock(d, n, IMAG, seed, psi_in)
    if nargin < 1, d = 64; end
    if nargin < 2, n = 2; end
    if nargin < 3, IMAG = 1; end
    % seed: optional; omitted or [] leaves the RNG untouched (random psi)
    if nargin >= 4 && ~isempty(seed), rng(seed); end

    % Paths relative to this file, so it works from any working directory
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, 'util'));
    addpath(qclabDir());

    % Random normalized target state in C^(d^n), or R^(d^n) if IMAG = 0;
    % psi_in (optional) replaces it with a given state, normalized here
    if nargin >= 5 && ~isempty(psi_in)
      assert(numel(psi_in) == d^n, 'psi_in must have d^n entries');
      psi = psi_in(:) / norm(psi_in);
    else
      psi = randn(d^n, 1) + 1i * IMAG * randn(d^n, 1);
      psi = psi / norm(psi);
    end
    ORIGINAL_PSI = psi;

    tic;
    cs = makeClubSequence(d, n);

    % gates collected here and added to the circuits in one insert
    % (see circuitFromGates)
    gates = cell(1, numel(cs) + 1);   % reduction gates, in order
    k = 0;

    for term = cs
        [c, cv, t, V] = singleClubHouseholder(term, psi, d);
        VGate = qclab.qgates.MatrixGate(t-1, V);
        if c ~= -1
            VGate = qclab.qgates.ControlledGate(VGate, c-1, t-1, cv - '0');
        end
        psi = VGate.apply('R', 'N', n, psi, 0, d);
        % psi(abs(psi) < 1e-6) = 0
        k = k + 1;
        gates{k} = VGate;
    end

    % the reduced state is e^{i theta}|0>: the reduction ends by removing it
    globalPhase = psi(1,1);
    gates{end} = qclab.qgates.qudit.SubspaceGate( ...
      qclab.qgates.Phase(n-1, real(globalPhase), -imag(globalPhase)), [1, 0], n-1);
    % reduction circuit: psi -> |0>; preparation circuit: its adjoint, i.e.
    % the adjoints of the gates in reverse order (built directly:
    % QCircuit.ctranspose assigns its gates one by one, which costs O(K^2)
    % like push_back)
    circuit = circuitFromGates(n, d, gates);
    prep = circuitFromGates(n, d, cellfun(@ctranspose, fliplr(gates), ...
      'UniformOutput', false));
    b_time = toc;

    % |0...0> as a vector: the bitstring form uses base2dec, which fails for d > 36
    e0 = zeros(d^n, 1);
    e0(1) = 1;
    tic;
    res = prep.simulate(e0).states;
    s_time = toc;

    err = norm(res - ORIGINAL_PSI);

    if nargout == 0
        % Called as a script: show results instead of returning them
        fprintf('Bullock QR state preparation, n = %d qudits, d = %d\n', n, d);
        fprintf('  gates      : %d\n', circuit.nbObjects);
        fprintf('  error      : %.3e  (||res - psi||_2)\n', err);
        fprintf('  build time : %.4f s\n', b_time);
        fprintf('  sim time   : %.4f s\n', s_time);
        clear circuit err b_time s_time   % avoid printing "ans"
    end
end
