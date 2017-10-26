#!/usr/local/bin/julia5
include("Eirene0_3_5.jl")
using MultivariateStats.classical_mds, JLD

const SAMPLING = length(ARGS) > 0 ? ARGS[1] : "uniform"

mds(M;dim=2) = classical_mds(M,dim)

function torus!(r,R)
    rR = [r,R-r]
    mn = minimum(rR)
    a,b = (mn/6)./rR
    dt,dp = map(x->2*asin(x > 1 ? 1 : x),[a,b])
    ts,ps = map((d)->d/2:d:2pi,[dt,dp])
    M = [(t,p) for t=ts,p=ps]
    T = [[(R + r*cos(a[1])).*cos(a[2]) (R + r*cos(a[1])).*sin(a[2]) r*sin(a[1])] for a in M]
    [T[:][i][j] for j=1:3,i=1:length(T[:])]
end

function _torus!(r,R,n=500)
    seed = 2pi*rand(n,2)
    x = (R + r*cos.(seed[:,1])).*cos.(seed[:,2])
    y = (R + r*cos.(seed[:,1])).*sin.(seed[:,2])
    z = r*sin.(seed[:,1])
    [x y z].'
end

function dmat(P;dim=2,d=(x,y)->norm(x-y))
    n = dim > 0 ? size(P)[dim] : length(P)
    D = Array{Float32,2}(n,n)
    for i=1:n
        D[i,i] = 0
        for j=i+1:n
            dist = dim > 0 ? dim > 1 ? d(P[:,i],P[:,j]) : d(P[i,:],P[j,:]) : d(P[i],P[j])
            D[i,j] = D[j,i] = dist
        end
    end
    D
end

function testpts(d,pts,e)
    # nbhd(p) = abs(d - p) < e
    # res = nbhd(e) ? true : false
    # for pt in pts
    #     res = nbhd(pt) ? true : res
    # end
    # res
    rand(Bool)
end

function getreps(E,M,nfeatures=0)
    m,n = size(M)
    bcode = barcode(E)
    nbcodes = size(bcode)[1]
    nfeatures = nfeatures > 0 ? nfeatures : nbcodes
    L = [bcode[i,2] - bcode[i,1] for i=1:nbcodes]
    p = sortperm(L,rev=true)
    RP = zeros(m,n)
    edge_reps = []
    vertex_reps = []
    for c in p[1:nfeatures]
        push!(edge_reps,[])
        reps = classrep(E,class=c)
        push!(vertex_reps,sort(union([],reps)))
        for k=1:size(reps)[2]
            i,j = map(x->reps[x,k],(1,2))
            RP[i,j] = RP[j,i] = M[i,j]
            push!(edge_reps[end],(i,j))
        end
    end
    L,RP,edge_reps,vertex_reps,bcode[p,:]#map(x->RP[x,x],vertex_reps)
end

# PARAMS
_d = 3
_rc = sqrt(2*_d/(_d+1))/2
r,R = (0.3,0.9)
low,high = (0.0,_rc*2*(R+r))
mds_dim = 2

T = SAMPLING == "random" ? _torus!(r,R,500) : torus!(r,R)
n = size(T)[2]
println("\n$n points")
println("(r,R):\t($r,$R)\n")

# UNFILTERED
_M = _rc*dmat(T)
_E = eirene(_M,rowsare="distances",lowerlim=low)#,upperlim=high)
_L,_RP,_edge_reps,_vertex_reps,_bcode = getreps(_E,_M)
_P = mds(_RP,dim=mds_dim)

# PLOT INIT
import PyPlot
const plt = PyPlot
fig = plt.figure(1,figsize=(11,6))
plt.ion()

# INIT UNFILTERED
ax1 = plt.subplot(231);ax1[:plot]([0,1],[0,1],c="black",alpha=0.5)
ax2 = mds_dim > 2 ? plt.subplot(232,projection="3d") : plt.subplot(232)
ax2[:axis]("off"); ax2[:axis]("equal")
ax3 = plt.subplot(233,projection="3d");ax3[:axis]("off")
ax3[:set_xlim](-1,1); ax3[:set_ylim](-1,1); ax3[:set_zlim](-1,1)

ax4 = plt.subplot(234)
ax5 = mds_dim > 2 ? plt.subplot(235,projection="3d") : plt.subplot(235)
ax6 = plt.subplot(236,projection="3d")
plt.tight_layout()

# PLOT UNFILTERED
# ax1[:scatter]([0,0],[r,R-r],s=30,c="black",alpha=0.2)
ax1[:scatter](0,r,s=30,c="green",alpha=0.2,marker="+")
ax1[:scatter](0,R-r,s=30,c="red",alpha=0.2,marker="+")
ax1[:scatter](_bcode[:,1],_bcode[:,2],s=5,alpha=0.5)
ax1[:scatter](_bcode[2,1],_bcode[2,2],c="green",s=20,alpha=0.5)
ax1[:scatter](_bcode[1,1],_bcode[1,2],c="red",s=20,alpha=0.5)
if mds_dim > 2
    ax2[:scatter](_P[1,:],_P[2,:],_P[3,:],s=5,alpha=0.5)
else
    ax2[:scatter](_P[1,:],_P[2,:],s=5,alpha=0.5)
end
ax3[:scatter](T[1,:],T[2,:],T[3,:],s=5,alpha=0.5)
for i=1:2
    color = i == 1 ? "red" : "green"
    ax2[:scatter]([_P[k,_vertex_reps[i]] for k=1:mds_dim]...,s=20,c=color,alpha=0.5)
    for e in _edge_reps[i]
        X,Y,Z = map(x->[T[x,e[1]],T[x,e[2]]],(1,2,3))
        ax3[:plot](X,Y,Z,c=color,alpha=0.5)
    end
end

function run(PARAM=6)
    # FILTERED
    M = copy(_M)
    cntr = sum([T[:,i] for i=1:n])/n
    D = [norm(T[:,i] - cntr) for i=1:n]
    R_minus_r, R_plus_r = map(f->f(D),[minimum,maximum])
    pts = [(R_plus_r - R_minus_r)/2, R_minus_r]
    _e = _rc*maximum((minimum(pts)/PARAM)./pts)
    println("\npts:\t$(pts)")
    println("_e:\t$(_e)\n")

    edges = []
    for i=1:n
        for j=i+1:n
            if testpts(M[i,j],pts,_e)
                push!(edges,(i,j))
            else
                M[i,j] = M[j,i] = Inf
            end
        end
    end

    nedges = length(edges)
    nedges_all = Int(round((n^2-n)/2))
    pct_edges = round(1000*nedges/nedges_all)/10
    println("$(nedges) of $(nedges_all) edges (~$(pct_edges)%)")

    E = eirene(M,rowsare="distances",lowerlim=low)#,upperlim=high)
    L,RP,edge_reps,vertex_reps,bcode = getreps(E,_M)
    P = mds(RP,dim=mds_dim)

    bcode[:,2] = map(i->(bcode[i,2]==Inf ? bcode[i,1] : bcode[i,2]),1:size(bcode)[1])

    # INIT FILTERED
    map(ax->ax[:cla](),[ax4,ax5,ax6])
    ax4[:plot]([0,1],[0,1],c="black",alpha=0.5)
    ax5[:axis]("off"); ax5[:axis]("equal")
    ax6[:axis]("off");ax6[:set_xlim](-1,1); ax6[:set_ylim](-1,1); ax6[:set_zlim](-1,1)

    # PLOT FILTERED
    # ax4[:scatter]([0,0],[r,R-r],s=30,c="black",alpha=0.2)
    # ax4[:plot]([0,0],[r-_e/2,r+_e/2],c="green",alpha=0.2)
    ax4[:scatter](0,r,s=30,c="green",alpha=0.2,marker="+")
    # ax4[:plot]([0,0],[R-r-_e/2,R-r+_e/2],c="red",alpha=0.2)
    ax4[:scatter](0,R-r,s=30,c="red",alpha=0.2,marker="+")
    ax4[:scatter](bcode[:,1],bcode[:,2],s=5,alpha=0.5)
    ax4[:scatter](bcode[2,1],bcode[2,2],c="green",s=20,alpha=0.5)
    ax4[:scatter](bcode[1,1],bcode[1,2],c="red",s=20,alpha=0.5)
    # ax4[:scatter](zeros(size(bcode)[1]),L,s=5,alpha=0.5,marker="^")
    if mds_dim > 2
        ax5[:scatter](P[1,:],P[2,:],P[3,:],s=5,alpha=0.5)
    else
        ax5[:scatter](P[1,:],P[2,:],s=5,alpha=0.5)
    end
    ax6[:scatter](T[1,:],T[2,:],T[3,:],s=5,alpha=0.5)
    for i=1:2
        color = i == 1 ? "red" : "green"
        ax5[:scatter]([P[k,vertex_reps[i]] for k=1:mds_dim]...,s=20,c=color,alpha=0.5)
        for e in edge_reps[i]
            X,Y,Z = map(x->[T[x,e[1]],T[x,e[2]]],(1,2,3))
            ax6[:plot](X,Y,Z,c=color,alpha=0.5)
            # if mds_dim > 2
            #     _X,_Y,_Z = map(x->[P[x,e[1]],P[x,e[2]]],(1,2,3))
            #     ax5[:plot](_X,_Y,_Z,c=color,alpha=0.5)
            # else
            #     _X,_Y = map(x->[P[x,e[1]],P[x,e[2]]],(1,2))
            #     ax5[:plot](_X,_Y,c=color,alpha=0.5)
            # end
        end
    end
    # plt.savefig("tex/figures/$(n)$(SAMPLING)$(PARAM).pdf")
    plt.savefig("tex/figures/$(n)$(SAMPLING)half.pdf")

    Dict(:parameter=>PARAM,:critical=>pts,:e=>_e,:edges=>edges,:M=>M,
        :barcode=>bcode,:lifetimes=>L,:repmat=>RP,:mds=>P,
        :edge_reps=>edge_reps,:vertex_reps=>vertex_reps)
end

# DICTS=[run(p) for p=5:0.5:12]
# save("data/$(n)$(SAMPLING).jld", :D, DICTS)
# D = run(10)

D=run()
save("data/$(n)$(SAMPLING)half.jld", "D", D)
