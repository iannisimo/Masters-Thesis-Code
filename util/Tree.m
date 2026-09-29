classdef Tree
  properties
    n
    d
    tree
  end

  methods
    function obj = Tree(vec, d)
      N = numel(vec);
      % round: in floating point log(d^n)/log(d) can miss n (log2(243)/log2(3) < 5)
      obj.n = round(logd(N, d));
      assert(d^obj.n == N, 'Tree: numel(vec) must be a power of d');
      obj.d = d;
      obj.tree = cell(1, obj.n+1);
      obj.tree{1} = vec;
      for i = 2:obj.n+1
        vec_ = obj.tree{i-1};
        vec_ = reshape(vec_, d, []);
        k = size(vec_, 2);
        vec_ = mat2cell(vec_, d, ones(1, k));
        % norm of the children with the phase of the first one (1 if it
        % is zero, the same convention s = 1 as makeHouseholder)
        vec_ = cellfun(@(col) norm(col) * firstPhase(col(1)), vec_);
        obj.tree{i} = reshape(vec_, [], 1);
      end
    end

    function phase = globalPhase(obj)
      phase = obj.tree{end};
    end

    function [phi, ctrl, target, ctrlStates] = getVector(obj, depth, idx)
      phi = obj.tree{end-depth}((idx - 1)*obj.d + 1:idx * obj.d);
      ctrl = 0:depth-2;
      target = depth - 1;
      ctrlStates = dec2base_(idx - 1, obj.d, size(ctrl, 2));
    end

    function [reflectors, ctrl, target, ctrlsStates] = getReflectors(obj, depth)
      vec_ = obj.tree{end-depth};
      vec_ = reshape(vec_, obj.d, []);
      k = size(vec_, 2);
      vec_ = mat2cell(vec_, obj.d, ones(1, k));
      reflectors = cellfun(@makeHouseholder, vec_, 'UniformOutput', false);
      ctrl = 0:depth-2;
      target = depth - 1;
      ctrlsStates = {};
      if(~isempty(ctrl))
        ctrlsStates = cellfun(@(idx) dec2base_(idx, obj.d, size(ctrl, 2)), num2cell(0:k-1), 'UniformOutput', false);
      end
    end
  end

end

function p = firstPhase(z)
% phase z/|z| of the first child, 1 when z = 0
if z == 0
  p = 1;
else
  p = z / abs(z);
end
end
