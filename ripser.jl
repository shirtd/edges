# const _rc = sqrt(2*_d/(_d+1))/2

# lower = Dict{Symbol,Function}(:matrix=>(M)->[[M[i,j] for j=1:i-1] for i=1:size(M)[1]],
# 		:points=>(P)->[[_rc*norm(P[:,i] - P[:,j]) for j=1:i-1] for i=1:size(P)[2]])
lower(M) = [[M[i,j] for j=1:i-1] for i=1:size(M)[1]]
# ripsin = Dict(T=>(D)->replace("$(lower[T](D))","],[",",\n")[16:end-2] for T in [:points,:matrix])

# ripsin(M) = replace("$(lower(M))",r"Array{(*?),1}\[Float32\[|\],(*?)\[",",\n")[16:end-2]
ripsin(M) = replace(replace("$(lower(M))",r"(\],(.*?)\[)|(\]\])",",\n"),r"(Array\{(.*?),1\}\[(.*?)\[)","")
run!(e,s) = (print("\trunning ripser..."); println(e[2],s); close(e[2]); readstring(e[1]))
# ripser!(D,T=:point;dim=1) = @time run!(readandwrite(`./ripser --dim $dim`),ripsin[T](D))
ripser!(M;dim=1) = @time run!(readandwrite(`./ripser --dim $dim`),ripsin(M))

T_str = Dict(:points=>(D)->"$(size(D)[2]) points in R$(size(D)[1])",
		:matrix=>(D)->"$(size(D)[1]) point distance matrix")

function ripser(M;dim=1)
    println("\n> computing barcodes 0-$dim of $(size(M)[1]) point distance matrix")
    stdin = ripser!(M,dim=dim)
    n = parse(Int,match(r"(?<=distance matrix with )(.*)(?= points)",stdin).match)
    llim,ulim = map(s->parse(Float32,s),split(match(r"(?<=value range: \[)(.*)(?=\])",stdin).match,","))
    dims = map(s->parse(Int,s),matchall(r"(?<=persistence intervals in dim )(.*)(?=:)",stdin))
    parse_floats(strs) = map(x->(v=tryparse(Float64,x); isnull(v) ? Inf : get(v)),strs)
    str_to_pairs(str) = map(s->parse_floats(split(s,",")),matchall(r"(?<= \[)(.*)(?=\))",str))
    pairs = [str_to_pairs(str) for str in split(stdin, r"persistence intervals in dim .*?:")[2:end]]
    barcode = Dict(dims[j]=>[pair[i] for pair in pairs[j], i in 1:2] for j=1:length(dims))
	max = -Inf
	for dim in dims
		if length(barcode[dim]) > 0
			mx = maximum(filter(x->(x<Inf),barcode[dim]))
			max = mx > max ? mx : max
		end
	end
    Dict(:n=>n,:lims=>(llim,ulim),:max=>max,:dims=>dims,:barcode=>barcode,:M=>M)
end
