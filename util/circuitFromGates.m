function circuit = circuitFromGates(n, d, gates)
% CIRCUITFROMGATES  QCircuit on n qudits of dimension d with the gates of the
% cell array gates, in order, added with a single insert.
% QCircuit.push_back (and every indexed assignment into the circuit) copies
% and re-validates the whole gate list, so building a circuit gate by gate
% costs O(K^2) for K gates; one insert of the whole list costs O(K).
circuit = qclab.QCircuit(n, 0, d);
if ~isempty(gates)
  circuit.insert(1:numel(gates), [gates{:}]);
end
end
