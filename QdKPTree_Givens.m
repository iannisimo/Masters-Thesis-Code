function [circuit, err, b_time, s_time] = QdKPTree_Givens(d, n, IMAG, seed, psi_in)
  % Usable both ways:
  %   QdKPTree_Givens            % run as a "script" (Run button / matlab -batch QdKPTree_Givens), uses defaults
  %   [c, e, b, s] = QdKPTree_Givens(d, n, IMAG, seed, psi_in)
  % NOTE: the amplitudes are squared with power() and the signs folded into
  % the leaf rotations, so complex states (IMAG ~= 0) are not supported yet.
  if nargin < 1, d = 3; end
  if nargin < 2, n = 3; end
  if nargin < 3, IMAG = 0; end
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
  % generate the tree with the squares of psi on the leaves
  % and the sums of the children in the internal nodes;
  % keep the signs in a separate list for the leaves
  t = cell(1, n+1);
  s = cell(1, 1);

  % generate the d-ary tree from psi
  for i = 1:n
    t{i} = cell(1, d^(n-i));
    if (i == 1)
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

  t{n+1} = sqrt(sum(power(psi, 2))) / norm(psi);

  % convert into binary tree for givens
  for i = 1:n
    for j = 1:d^(n-i)
      dary = t{i}{j};
      dary = reshape(dary, 1, []);
      dary = [dary, zeros(1, (power(2, ceil(log2(d))) - numel(dary)))];
      sub_depth = ceil(log2(d));
      t{i}{j} = cell(1, sub_depth);
      for k = 1:sub_depth
        D = numel(dary);
        elems = D / (2^k);
        t{i}{j}{k} = cell(1, elems);
        for l = 1:elems
          if k == 1
            t{i}{j}{k}{l} = dary((l-1)*2 + 1:l*2);
          else
            for m = 1:2
              t{i}{j}{k}{l}(m) = sum(t{i}{j}{k-1}{(l-1) * 2 + m});
            end
          end
        end
      end
    end
  end

  % generate circuit from QdKP-Tree; gates collected here and added to the
  % circuit in one insert (see circuitFromGates)
  gates = {};

  for i = n:-1:1
    for j = 1:numel(t{i})
      % set controls and target
      if i < n
        ctrl = 0:n-i-1;
        ctrlStates = dec2base_(j-1, d, n-i);
      end
      target = n-i;
      sub_depth = numel(t{i}{j});
      for k = sub_depth:-1:1
        for l = 1:numel(t{i}{j}{k})
          old = 1;
          if (k < sub_depth)
            old = t{i}{j}{k+1}{ceil(l/2)}(mod(l-1, 2) + 1);
          elseif i < n
            child_idx = mod(j-1, d) + 1;
            old = t{i+1}{ceil(j/d)}{1}{ceil(child_idx/2)}(mod(child_idx-1, 2) + 1);
          end
          % givens subspace calculation
          leaf_l = (l-1) * 2;
          leaf_r = leaf_l + 1;
          sub_l = leaf_l * 2^(k-1);
          sub_r = leaf_r * 2^(k-1);
          phi = t{i}{j}{k}{l};
          cs = sqrt(phi(1));
          sn = sqrt(phi(2));
          % leaf level: fold the signs of psi into the rotation
          if i == 1 && k == 1 && sub_r < d
            cs = real(s{1}{j}(sub_l+1)) * cs;
            sn = real(s{1}{j}(sub_r+1)) * sn;
          elseif i == 1 && sub_r == d-1
            % odd d: leaf d-1 has a padded partner, so its sign goes on the
            % first rotation where its subtree (padded apart from it) is
            % on the right
            sn = real(s{1}{j}(d)) * sn;
          end
          % skip only when the rotation is the identity
          if sn == 0 && cs >= 0, continue; end
          ry = qclab.qgates.RotationY(target, cs / sqrt(old), sn / sqrt(old), true);
          % SubspaceGate embeds the qubit gate in a d x d unitary
          % and places the original operator on the subspace defined
          % as the second parameter
          ry = qclab.qgates.qudit.SubspaceGate(ry, [sub_l, sub_r], target);
          if i < n
            ry = qclab.qgates.MControlledGate(ry, ctrl, target, ctrlStates);
          end
          gates{end+1} = ry; %#ok<AGROW>
        end
      end
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
    fprintf('QdKP-Tree (Givens) state preparation, n = %d qudits, d = %d, IMAG = %d\n', n, d, IMAG);
    fprintf('  gates      : %d\n', circuit.nbObjects);
    fprintf('  error      : %.3e  (||res - psi||_2)\n', err);
    fprintf('  build time : %.4f s\n', b_time);
    fprintf('  sim time   : %.4f s\n', s_time);
    clear circuit err b_time s_time   % avoid printing "ans"
  end
end
