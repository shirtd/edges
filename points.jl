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
    nbhd(p) = abs(d - p) < e
    res = nbhd(e) ? true : false
    for pt in pts
        res = nbhd(pt) ? true : res
    end
    res
end

function testwindow(d,w)
    w[1] <= d && d <= w[2]
end
