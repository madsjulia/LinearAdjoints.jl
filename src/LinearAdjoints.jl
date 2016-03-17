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

end
