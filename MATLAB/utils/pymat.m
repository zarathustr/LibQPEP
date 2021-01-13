function out = pymat(obj, m, n)
out = double(py.array.array('d', py.numpy.nditer(obj))); %d is for double, see link below on types
out = reshape(out, [m n]);
end