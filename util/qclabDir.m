function dir = qclabDir()
% QCLABDIR  Folder of the QCLAB fork (github.com/iannisimo/qclab, branch
% qudit): the qclab submodule at the root of this repo.
  dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'qclab');
  if ~exist(fullfile(dir, '+qclab'), 'dir')
    error('qclabDir:notFound', ['QCLAB not found in %s; fetch the submodule:\n' ...
          '  git submodule update --init'], dir);
  end
end
