function result = ndarray2mat(ndarray)
% Convert an ndarray from numpy to a MATLAB array
data_size = cellfun(@int64,cell(ndarray.shape));
if length(data_size) == 1
    % Scalars are trivial
    result = double(py.array.array('d', py.numpy.nditer(ndarray)));
elseif length(data_size) == 2
    % order='F' is used to get data in column-major order (as in Fortran
    % 'F' and MATLAB)
    result = reshape(double(py.array.array('d', ...
        py.numpy.nditer(ndarray, pyargs('order', 'F')))), ...
        data_size);
else
    % For multidimensional arrays more manipulation is required
    % First recover in Python order (C contiguous order)
    result = double(py.array.array('d', ...
        py.numpy.nditer(ndarray, pyargs('order', 'C'))));
    % Switch the order of the dimensions (as Python views this in the
    % opposite order to MATLAB) and reshape to the corresponding C-like
    % array
    result = reshape(result,fliplr(data_size));
    % Now transpose rows and columns of the 2D sub-arrays to arrive at the
    % correct MATLAB structuring
    result = permute(result, [length(data_size):-1:1]);
end
