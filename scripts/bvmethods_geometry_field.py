"""Private geometry-only field diagnostic; does not select a solver profile.

Coordinates set the inverse-length unit of the returned dE/dpsi operator.
Box transfer is limited to unweighted manifold Delaunay Tri3 meshes with
non-obtuse material boundaries. Shared-edge weights include both materials;
the field stencil uses only neighbors of the requested material region.
"""
import collections
import math
import numpy as np


def delaunay_boxes(points, triangles, regions):
    p=np.asarray(points,float)
    if p.ndim!=2 or p.shape[1]!=2 or not np.all(np.isfinite(p)):
        raise ValueError('finite two-dimensional coordinates required')
    if len(triangles)!=len(regions) or not len(triangles):
        raise ValueError('matching nonempty triangle regions required')
    edges=collections.defaultdict(list);kappa=np.zeros((len(triangles),3))
    for c,ns in enumerate(triangles):
        if len(ns)!=3 or len(set(ns))!=3 or min(ns)<0 or max(ns)>=len(p):
            raise ValueError('three distinct valid vertices required')
        for k in range(3):
            i,j,o=ns[k],ns[(k+1)%3],ns[(k+2)%3]
            # Extended arithmetic is used for geometry, not fitted to native E.
            u=np.asarray(p[i],np.longdouble)-p[o];v=np.asarray(p[j],np.longdouble)-p[o]
            cross=abs(u[0]*v[1]-u[1]*v[0])
            if not cross>0:raise ValueError('degenerate triangle')
            value=float(u@v/(2*cross))
            if not math.isfinite(value):raise ValueError('nonfinite cotangent')
            edges[tuple(sorted((i,j)))].append((c,k,value))
    transfers=0
    for parts in edges.values():
        if len(parts)>2:raise ValueError('nonmanifold edge')
        total=math.fsum(v for c,k,v in parts);tol=128*np.finfo(float).eps*max(1.,*(abs(v) for c,k,v in parts))
        if total < -tol:raise ValueError('globally non-Delaunay edge')
        negative=min(parts,key=lambda x:x[2])
        if negative[2]<0:
            same=len(parts)==2 and regions[parts[0][0]]==regions[parts[1][0]]
            if not same and negative[2]<-tol:raise ValueError('obtuse material boundary')
            for c,k,v in parts:kappa[c,k]=0 if (c,k)==negative[:2] else max(total,0.) if same else v
            transfers+=int(same)
        else:
            for c,k,v in parts:kappa[c,k]=v
    volumes=np.zeros_like(kappa)
    for c,ns in enumerate(triangles):
        products=[kappa[c,k]*np.sum((p[ns[(k+1)%3]]-p[ns[k]])**2) for k in range(3)]
        for k in range(3):volumes[c,k]=(products[k]+products[(k-1)%3])/4
    return kappa,volumes,transfers


def region_field_stencil(points, triangles, regions, target_region):
    """Generate a region stencil, or all mesh neighbors for target_region=None."""
    p=np.asarray(points,float);_,volumes,transfers=delaunay_boxes(p,triangles,regions)
    weights=collections.defaultdict(float);adj=collections.defaultdict(set)
    for c,ns in enumerate(triangles):
        for k,i in enumerate(ns):
            for l,j in enumerate(ns):
                if i==j:continue
                weights[i,j]+=(volumes[c,k]+volumes[c,l])/np.sum((p[j]-p[i])**2)
                if target_region is None or regions[c]==target_region:adj[i].add(j)
    if not adj:raise ValueError('empty target region')
    entries=[]
    for i,neighbors in sorted(adj.items()):
        ns=sorted(neighbors);d=p[ns]-p[i];w=np.array([weights[i,j] for j in ns]);mat=(d.T*w)@d
        if np.linalg.matrix_rank(mat)<2:raise ValueError('rank-deficient field neighborhood')
        A=-np.linalg.solve(mat,d.T*w)
        for k,j in enumerate(ns):entries.append((i,j,A[:,k]))
        entries.append((i,i,-np.sum(A,axis=1)))
    return entries,dict(transferred_edges=transfers,region_nodes=len(adj))


def semiconductor_node_field_stencil(points, triangles, regions, semiconductor_region):
    """BVmethods candidate: global recovery, overridden at Si-region vertices.

    This single-semiconductor diagnostic does not qualify heterojunctions.
    Shared insulator datasets inherit the same field at a Si vertex.
    """
    global_entries,identity=region_field_stencil(points,triangles,regions,None)
    si_entries,si_identity=region_field_stencil(points,triangles,regions,semiconductor_region)
    nodes={i for i,j,a in si_entries}
    return [e for e in global_entries if e[0] not in nodes]+si_entries,dict(
        transferred_edges=identity['transferred_edges'],global_nodes=identity['region_nodes'],
        semiconductor_nodes=si_identity['region_nodes'])
