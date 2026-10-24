function result = mat2ndarray(mat)
% Convert a MATLAB array into an ndarray
data_size = size(mat);
if length(data_size) == 1
    % Scalars are trivial
    result=py.numpy.array(mat);
elseif length(data_size) == 2
    % A transpose operation is required either in MATLAB, or in Python due
    % to the difference between row major and column major ordering
    transpose = mat';
    % Pass the array to Python as a vector, and then reshape to the correct
    % size
    result = py.numpy.reshape(transpose(:)', int32(data_size));
    % Flatten out singleton dimension
    if any(data_size == 1)
        result = result.flatten();
    end
else
    % For an n-dimensional array, transpose the first two dimensions to
    % sort the storage ordering issue
    transpose = permute(mat, [length(data_size):-1:1]);
    % Pass it to Python, and then reshape to the Python style of matrix
    % sizing
    result = py.numpy.reshape(transpose(:)', int32(fliplr(size(transpose))));
end
