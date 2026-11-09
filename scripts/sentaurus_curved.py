"""Certified straight-chord representation of one convex Rectangle fillet."""
import math

CHORD_ERROR_UM = 1e-5


def polygon(body):
    geometry = body.get('geometry', body)
    x0,y0 = geometry['lower_left'];x1,y1 = geometry['upper_right']
    corners = [(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
    fillets = geometry.get('fillets', [])
    points=[]
    for i,corner in enumerate(corners):
        selected=[f for f in fillets if tuple(f['corner'])==corner]
        if not selected:
            points.append(corner);continue
        f=selected[0];r=f['radius_um']
        previous,next_=corners[(i-1)%4],corners[(i+1)%4]
        def toward(p):
            length=math.dist(p,corner)
            return ((p[0]-corner[0])/length,(p[1]-corner[1])/length)
        u,v=toward(previous),toward(next_)
        center=(corner[0]+r*(u[0]+v[0]),corner[1]+r*(u[1]+v[1]))
        start=(corner[0]+r*u[0],corner[1]+r*u[1])
        a=math.atan2(start[1]-center[1],start[0]-center[0])
        n=max(1,math.ceil((math.pi/2)/(2*math.acos(1-min(CHORD_ERROR_UM/r,1.)))))
        points.extend((center[0]+r*math.cos(a+math.pi*j/(2*n)),center[1]+r*math.sin(a+math.pi*j/(2*n))) for j in range(n+1))
    # A radius may consume a complete short side; merge coincident tangent vertices.
    clean=[]
    for point in points:
        if not clean or math.dist(point,clean[-1])>1e-12:clean.append(point)
    if len(clean)>1 and math.dist(clean[0],clean[-1])<1e-12:clean.pop()
    return clean


def area(body):
    p=polygon(body)
    return sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(p,p[1:]+p[:1]))/2


def contains(body,point,tol=1e-9):
    p=polygon(body)
    return all((b[0]-a[0])*(point[1]-a[1])-(b[1]-a[1])*(point[0]-a[0]) >= -tol*math.dist(a,b)
               for a,b in zip(p,p[1:]+p[:1]))
