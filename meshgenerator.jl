using Plots, PyCall
using SparseArrays
using TriplotRecipes
function delaunay(p)
    tri = pyimport("matplotlib.tri")
    t = tri[:Triangulation](p[:,1], p[:,2])
    return Int64.(t[:triangles] .+ 1)
end
using LinearAlgebra

edgetol=1e-7
A=0.4
T=0.5
epsilon=1e-12
function all_edges(t)
    etag = vcat(t[:,[1,2]], t[:,[2,3]], t[:,[3,1]])
    etag = hcat(sort(etag, dims=2), 1:3*size(t,1))
    etag = sortslices(etag, dims=1)
    dup = all(etag[2:end,1:2] - etag[1:end-1,1:2] .== 0, dims=2)[:]
    keep = .![false;dup]
    edges = etag[keep,1:2]
    emap = cumsum(keep)
    invpermute!(emap, etag[:,3])
    emap = reshape(emap,:,3)
    dup = [dup;false]
    dup = dup[keep]
    bndix = findall(.!dup)
    return edges, bndix, emap
end


"""
    e = boundary_nodes(t)

Find all boundary nodes in the triangulation `t`.
"""
function boundary_nodes(t)
    edges, boundary_indices, _ = all_edges(t)
    return unique(edges[boundary_indices,:][:])
end

"""
    tplot(p, t, u=nothing)

If `u` == nothing: Plot triangular mesh with nodes `p` and triangles `t`.
If `u` == solution vector: Plot filled contour color plot of solution `u`.
"""
function tplot(p, t, u=nothing)
    if u == nothing
        p = trimesh(p[:,1], p[:,2], t', linecolor=:black,
            aspect_ratio=:equal, color=RGB(0.8,1,0.8))
    else
        p = tripcolor(p[:,1], p[:,2], u, t', color=:viridis,
            aspect_ratio=:equal)
    end
    p
end

"""
    inside = inpolygon(p, pv)

Determine if each point in the N x 2 node array `p` is inside the polygon
described by the NE x 2 node array `pv`.
"""
function inpolygon(p, pv)
    if ndims(p) == 2 && size(p,2) == 2
        return [ inpolygon(p[ip,:], pv) for ip = 1:size(p,1) ]
    end
    cn = 0
    for i = 1:size(pv,1) - 1
        if pv[i,2] <= p[2] && pv[i+1,2] > p[2] ||
           pv[i,2] > p[2] && pv[i+1,2] <= p[2]
            vt = (p[2] - pv[i,2]) / (pv[i+1,2] - pv[i,2])
            if p[1] < pv[i,1] + vt * (pv[i+1,1] - pv[i,1])
                cn += 1
            end
        end
    end
    return cn % 2 == 1
end


function pointgenerator(pv,hmax) #this function is supposed to do part b for us that is create nodes along each line segment 
    p=[-Inf -Inf]
    occured=[[-Inf,-Inf]] 
    countofvertices=0
    for i in 1:size(pv,1)-1
        dist=norm(pv[i+1,:]-pv[i,:])
        for j in 0.0:hmax:dist
            if !(pv[i,:]+j/dist*(pv[i+1,:]-pv[i,:]) in occured)
                p=[p; [pv[i,1]+j/dist*(pv[i+1,1]-pv[i,1]) pv[i,2]+j/dist*(pv[i+1,2]-pv[i,2])]]
                push!(occured, pv[i,:]+j/dist*(pv[i+1,:]-pv[i,:]))
                countofvertices+=1
            end
        end
            
    end

p=p[2:end,:]
distancelist=[norm(p[1+mod(i,size(p,1)),:]-p[i,:]) for i=1:size(p,1)]
p=[[p[i,1] for i in 1:size(p,1) if distancelist[i]>=edgetol] [p[i,2] for i in 1:size(p,1) if distancelist[i]>=edgetol]]
distancelist=[norm(p[1+mod(i,size(p,1)),:]-p[i,:]) for i=1:size(p,1)]
return p
end



function triArea(tri) #does the getting a triangle area. 
return abs((tri[2,1]+tri[1,1])*(tri[2,2]-tri[1,2])+(tri[3,1]+tri[2,1])*(tri[3,2]-tri[2,2])+(tri[1,1]+tri[3,1])*(tri[1,2]-tri[3,2]))/2
end

function tri_circumcenter(tri) #gets circumcenter of each triangle
X=tri[:,1]
Y=tri[:,2]
V=[1,1,1]
L=@.X^2+Y^2
M_y=det([V L Y])
M_x=det([V X L])
M_d=det([V X Y])
return ([M_y,M_x]/(2*M_d))
end
function tri_areas(p, t)
tri_pts = p[t, :]
 return 0.5 * abs.(tri_pts[:, 1, 1] .*
 (tri_pts[:, 2, 2] - tri_pts[:, 3, 2]) +
 tri_pts[:, 2, 1] .* (tri_pts[:, 3, 2] - tri_pts[:, 1, 2]) +
tri_pts[:, 3, 1] .* (tri_pts[:, 1, 2] - tri_pts[:, 2, 2]))
end
indicestopoints(p,l)=[p[l[1],1] p[l[1],2];p[l[2],1] p[l[2],2]; p[l[3],1 ] p[l[3],2]] #converts the list of indices into a list of Points
tri_centroid(tri)=[sum(tri[:,1])/3,sum(tri[:,2])/3]
function triAreas(p,t)
tri_pts=p[t,:]
    return [triArea(tri_pts[i,:,:]) for i=1:size(tri_pts,1)]

end
function triangle_generator_and_remover(pv,p_orig) #this function does steps c and d for us in our code
    triangulationobj=delaunay(pv)
    t=[0 0 0 ]
    k=size(triangulationobj,1)
    for i in 1:k
        centroid=(tri_centroid(indicestopoints(pv,triangulationobj[i,:])))
        area=triArea(indicestopoints(pv,triangulationobj[i,:]))
        if (inpolygon(centroid,p_orig) && area>=epsilon)
            t=[t; [triangulationobj[i,1] triangulationobj[i,2] triangulationobj[i,3]]]
        end
    end

    return t[2:end,:]
end

function pmesh(pv,hmax,nref=0)
    p=pointgenerator(pv,hmax)
    t=triangle_generator_and_remover(p,pv)
    arealist=[triArea(indicestopoints(p,t[i,:])) for i in 1:size(t,1)]
    maxarea=maximum(arealist)
    coordofmaxarea=findfirst(arealist.==maxarea)
    if (maxarea>hmax^2/2)
        newpointtopass=tri_circumcenter(indicestopoints(p,t[coordofmaxarea,:]))
        p=[p; [newpointtopass[1] newpointtopass[2]]]
    end
    while (maxarea>hmax^2/2)
        t=triangle_generator_and_remover(p,pv)
        arealist=[triArea(indicestopoints(p,t[i,:])) for i in 1:size(t,1)]
        maxarea=maximum(arealist)
        coordofmaxarea=findfirst(arealist.==maxarea)
        if (maxarea>hmax^2/2)
            newpointtopass=tri_circumcenter(indicestopoints(p,t[coordofmaxarea,:]))
            p=[p; [newpointtopass[1] newpointtopass[2]]]
        end
    end
    


for i in 1:nref
edges,boundary_nodes,_=all_edges(t)
midpoints=dropdims(sum(p[edges,:],dims=2)/2,dims=2)
p=[p;midpoints]
t=triangle_generator_and_remover(p,pv)

end
e=boundary_nodes(t)
return p,t,e

end

function subdivide(p,pv,t)
edges,boundary_nodes,_=all_edges(t)
midpoints=dropdims(sum(p[edges,:],dims=2)/2,dims=2)
p=[p;midpoints]
t=triangle_generator_and_remover(p,pv)
    return p,t
end
