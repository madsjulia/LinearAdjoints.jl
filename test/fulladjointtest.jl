import FDDerivatives
@LinearAdjoints.assemblevector (k, f) b function rhs(k, f)
	n = length(f)
	b = Array(Float64, length(f))
	for i = 1:length(f)
		b[i] = f[i] / n
	end
	return b
end
@LinearAdjoints.assemblesparsematrix (k, f) x function laplacian(k, f)
	n = length(f)
	I = Int[]
	J = Int[]
	V = Float64[]
	dx = 1 / (n - 1)
	for i = 1:n - 1
		LinearAdjoints.addentry(I, J, V, i, i, -2 * k / dx)
		LinearAdjoints.addentry(I, J, V, i + 1, i, k / dx)
		LinearAdjoints.addentry(I, J, V, i, i + 1, k / dx)
	end
	LinearAdjoints.addentry(I, J, V, n, n, -2 * k / dx)
	return sparse(I, J, V)
end
const k0 = float(pi)
xs = linspace(0, 1, length(h) + 2)[2:end - 1]
const hobs = randn(length(xs))
function objfunc(x)
	k = x[1]
	f = x[2:end]
	h, _ = handgrad(k, f)
	return objfunc(h, k, f)
end
objfunc_fdgrad = FDDerivatives.makegradient(objfunc)
function objfunc(h, k, f)
	return sum((h - hobs) .^ 2) + (k - k0) ^ 2
end
function objfunc_h(h, k, f)
	return 2 * (h - hobs)
end
function objfunc_p(h, k, f)
	result = zeros(1 + length(f))
	result[1] = 2 * (k - k0)
	return result
end
@LinearAdjoints.adjoint handgrad laplacian rhs objfunc_h objfunc_p
n = 10
k = k0
gwsink = fill(-2 * k0, n)
h, grad = handgrad(k, gwsink)
fdgrad = objfunc_fdgrad([k; gwsink])
@test_approx_eq_eps fdgrad grad sqrt(eps(Float64)) * 100
