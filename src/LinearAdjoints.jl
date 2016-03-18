module LinearAdjoints

import MetaProgTools

const specialsymbolI = :___la___I
const specialsymbolJ = :___la___J
const specialsymbolV = :___la___V

include("transforms.jl")
include("calculus.jl")
include("sparsematrix.jl")
include("vector.jl")

function writegetlinearindex(vars::Vector)
	gli = quote end
	for i = 1:length(vars)
		priorlengthexpression = :(0)
		for j = 1:i - 1
			priorlengthexpression = :($priorlengthexpression + length($(vars[j])))
		end
		q = quote
			function getlinearindex(::Type{$(Val{vars[i]})})
				return $priorlengthexpression + 1
			end
			function getlinearindex(::Type{$(Val{vars[i]})}, indices...)
				linearindex = $priorlengthexpression + 1
				for i = 1:length(indices)
					offset = 1
					for j = 1:i - 1
						offset *= size($(vars[i]), j)
					end
					linearindex += (indices[i] - 1) * offset
				end
				return linearindex
			end
		end
		push!(gli.args, q)
	end
	return gli
end

macro adjoint(name, assemble_A_func_symbol, assemble_b_func_symbol, objfunc_x_symbol, objfunc_p_symbol)
	q = quote
		function $name(args...)
			A = $assemble_A_func_symbol(args...)
			b = $assemble_b_func_symbol(args...)
			A2 = copy(A)
			Af = factorize(A)
			x = Af \ b
			g_x = $objfunc_x_symbol(x, args...)
			lambda = Af' \ g_x
			A_px = $(parse(string(assemble_A_func_symbol, "_px")))(x, args...)
			b_p = $(parse(string(assemble_b_func_symbol, "_p")))(args...)
			gradient = (b_p - A_px) * lambda + $objfunc_p_symbol(x, args...)
			return x, gradient
		end
	end
	return :($(esc(q)))
end

end
