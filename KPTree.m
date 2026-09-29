function [circuit, err, b_time, s_time] = KPTree(d, n, IMAG, seed, psi_in)
  % Usable both ways:
  %   KPTree            % run as a "script" (Run button / matlab -batch KPTree), uses defaults
  %   [c, e, b, s] = KPTree(d, n, IMAG, seed, psi_in)
  % NOTE: qubits only (d = 2), and the amplitudes are squared with power(), not
  % abs()^2, so complex states (IMAG ~= 0) are not supported yet.
  if nargin < 1, d = 2; end
  if nargin < 2, n = 2; end
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
  % generate the tree with the square of the elements
  % of psi on the leaves, and the sum of the children
  % on the internal nodes; keep sign in a separate ordered
  % list corresponding to the leaves of t
  t = cell(1, n+1);
  s = cell(1,1);

  for i = 1:n
    t{i} = cell(1, d^(n-i));
    if(i == 1)
      s{i} = cell(1, d^(n-i));
    end
    for j = 1:d^(n-i)
      if i == 1
        start = (j-1)*d;
        end_ = start + d - 1;
        t{i}{j} = power(psi(start+1:end_+1),2);
        s{i}{j} = sign(real(psi(start+1:end_+1))) + 1i * sign(imag(psi(start+1:end_+1)));
      else
        for k = 1:d
          t{i}{j}(k) = sum(t{i-1}{(j-1)*d + k});
        end
      end
    end
  end

  t{n+1} = psi(1) / norm(psi(1));

  % gates collected here and added to the circuit in one insert
  % (see circuitFromGates)
  gates = {};

  % phase = t{n+1};
  % circuit.push_back(qclab.qgates.PauliX(0));
  % circuit.push_back(qclab.qgates.Phase(0, real(phase), imag(phase)));
  % circuit.push_back(qclab.qgates.PauliX(0));

  % foreach pair in the tree, generate the
  % corresponding ry. For the last step,
  % add the sign data
  for i = n:-1:1
    for j = 1:numel(t{i})
      old = 1;
      if i < n
        old = t{i+1}{ceil(j/2)}(mod(j-1, 2) + 1);
      end
      phi = t{i}{j};
      % a node with no probability only acts on zero amplitudes: skip it
      % (its rotation would be 0/0)
      if old == 0, continue; end
      cs = 1;
      ss = 1;
      if i == 1
        cs = s{i}{j}(1);
        ss = s{i}{j}(2);
      end
      ry = qclab.qgates.RotationY(n-i, cs*sqrt(phi(1))/sqrt(old), ss*sqrt(phi(2))/sqrt(old), true);
      if i < n
        ctrl = 0:n-i-1;
        target = n-i;
        ctrlStates = dec2base_(j-1, d, n-i);
        ry = qclab.qgates.MControlledGate(ry, ctrl, target, ctrlStates);
      end
      gates{end+1} = ry; %#ok<AGROW>
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
    fprintf('KP-Tree state preparation, n = %d qudits, d = %d, IMAG = %d\n', n, d, IMAG);
    fprintf('  gates      : %d\n', circuit.nbObjects);
    fprintf('  error      : %.3e  (||res - psi||_2)\n', err);
    fprintf('  build time : %.4f s\n', b_time);
    fprintf('  sim time   : %.4f s\n', s_time);
    clear circuit err b_time s_time   % avoid printing "ans"
  end
end
