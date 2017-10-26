#!/usr/local/bin/julia5 -i
const _d = 3
const _rc = sqrt(2*_d/(_d+1))/2

include("points.jl")
include("ripser.jl")
include("plot.jl")

r,R = (0.3,0.9)
low,high = (0.0,_rc*2*(R+r))

# T = torus!(r,R)
T = _torus!(r,R,1000)
n = size(T)[2]

nedges_all = Int(round((n^2-n)/2))

println("\n$n points")
println("(r,R):\t($r,$R)\n")

# UNFILTERED
_M = _rc*dmat(T)
_R = ripser(_M)
plotR(ax1,_R)

ax2[:plot]([0,1.2*_R[:max]],[0,1.2*_R[:max]],c="black",alpha=0.5)
# ax2[:plot]([0,1.0],[0,1.0],c="black",alpha=0.5)

# # FILTERED
# M = copy(_M)
# cntr = sum([T[:,i] for i=1:n])/n
# D = [norm(T[:,i] - cntr) for i=1:n]
# R_minus_r, R_plus_r = map(f->f(D),[minimum,maximum])
# pts = [(R_plus_r - R_minus_r)/2, R_minus_r]
# _e = _rc*maximum((minimum(pts)/4)./pts)
#
# println("\npts:\t$(pts)")
# println("_e:\t$(_e)")
#
# edges = []
# for i=1:n
#     for j=i+1:n
#         if rand(Bool)#testpts(M[i,j],pts,_e)
#             push!(edges,(i,j))
#         else
#             M[i,j] = M[j,i] = Inf
#         end
#     end
# end
#
# nedges = length(edges)
# pct_edges = round(1000*nedges/nedges_all)/10
# println("$(nedges) of $(nedges_all) edges (~$(pct_edges)%)")
#
# R = ripser(M)
# plotR(ax2,R)

frames = 10
overlap = 0.5
_e = high/frames/overlap
windows = [(_e*(i*overlap-1),_e*i*overlap) for i=1:frames]

println("$frames windows")
println("overlap: $overlap")
println("_e: $_e")

edges = [[] for window in windows]
Ms = [copy(_M) for window in windows]
for i=1:n
    for j=i+1:n
        d = _M[i,j]
        for k=1:frames
            if d < _e || testwindow(d,windows[k])
                push!(edges[k],[i,j])
            else
                Ms[k][i,j] = Ms[k][j,i] = Inf
            end
        end
    end
end


# Rs = [ripser(M) for M in Ms]

# for k=1:frames
#     nedges = length(edges[k])
#     pct_edges = round(1000*nedges/nedges_all)/10
#     println("window $k: $(nedges) of $(nedges_all) edges (~$(pct_edges)%)")
#     R = ripser(Ms[k])
#     plotR(ax2,R,1,reset=false)
#     print("...")
#     readline(STDIN)
# end
