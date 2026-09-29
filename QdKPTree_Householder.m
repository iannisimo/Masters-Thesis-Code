function [circuit, err, b_time, s_time] = QdKPTree_Householder(d, n, IMAG, seed, psi_in)
  % Usable both ways:
  %   QdKPTree_Householder            % run as a "script" (Run button / matlab -batch QdKPTree_Householder), uses defaults
  %   [c, e, b, s] = QdKPTree_Householder(d, n, IMAG, seed, psi_in)
  if nargin < 1, d = 3; end
  if nargin < 2, n = 3; end
  if nargin < 3, IMAG = 1; end
  % seed: optional; omitted or [] leaves the RNG untouched (random psi)
  if nargin >= 4 && ~isempty(seed), rng(seed); end

  % Paths relative to this file, so it works from any working directory
  here = fileparts(mfilename('fullpath'));
  addpath(fullfile(here, 'util'));
  addpath(qclabDir());

  % generate random normalized psi, in C^(d^n) or in R^(d^n) if IMAG = 0;
  % psi_in (optional) replaces it with a given state, normalized here
  if nargin >= 5 && ~isempty(psi_in)
    assert(numel(psi_in) == d^n, 'psi_in must have d^n entries');
    psi = psi_in(:) / norm(psi_in);
  else
    psi = randn(d^n, 1) + 1i * IMAG * randn(d^n, 1);
    psi = psi / norm(psi);
  end

  tic;
  % generate the tree with psi on the leaves
  % and the (phased) norms of the children in the
  % internal nodes
  t = Tree(psi, d);

  % gates collected here and added to the circuit in one insert
  % (see circuitFromGates)

  % add initial phase gate with element from root
  phase = t.globalPhase();
  phase = qclab.qgates.Phase(n-1, real(phase), imag(phase));
  phase = qclab.qgates.qudit.SubspaceGate(phase, [1, 0], n-1);
  gates = {phase};

  % foreach |v> in the tree, generate the householder
  % and apply it in the circuit, with its corresponding
  % controls
  for i = 1:n
    [ref, ctrl, target, ctrlsStates] = t.getReflectors(i);
    for j = 1:numel(ref)
      HRGate = qclab.qgates.MatrixGate(target, ref{j}');
      if ~isempty(ctrl)
        HRGate = qclab.qgates.MControlledGate(HRGate, ctrl, target, ctrlsStates{j});
      end
      gates{end+1} = HRGate; %#ok<AGROW>
    end
  end
  circuit = circuitFromGates(n, d, gates);
  b_time = toc;

  % |0...0> as a vector: the bitstring form uses base2dec, which fails for d > 36
  e0 = zeros(d^n, 1);
  e0(1) = 1;
  tic;
  res = circuit.simulate(e0).states;
  s_time = toc;

  % test whether the state of the circuit corresponds
  % to the original psi
  err = norm(res - psi);

  if nargout == 0
    % Called as a script: show results instead of returning them
    fprintf('QdKP-Tree (Householder) state preparation, n = %d qudits, d = %d, IMAG = %d\n', n, d, IMAG);
    fprintf('  gates      : %d\n', circuit.nbObjects);
    fprintf('  error      : %.3e  (||res - psi||_2)\n', err);
    fprintf('  build time : %.4f s\n', b_time);
    fprintf('  sim time   : %.4f s\n', s_time);
    clear circuit err b_time s_time   % avoid printing "ans"
  end
end
