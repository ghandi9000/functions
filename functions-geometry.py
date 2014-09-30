## Solid angle subtended by triangle
from scipy import dot, arctan2, pi
from scipy.linalg import norm, det


def tri_projection(a, b, c):
    """Given three 3-vectors, a, b, and c."""
    determ = det((a, b, c))
    
    al = norm(a)
    bl = norm(b)
    cl = norm(c)
  
    div = al*bl*cl + dot(a, b)*cl + dot(a, c)*bl + dot(b, c)*al
    at = arctan2(determ, div)
    if at < 0:
        at += pi  # If det > 0 and div < 0 arctan2 returns < 0, so add pi.
    omega = 2 * at
  
    return omega

    
## Testing
a = (1, 0, 0)
b = (0, 1, 0)
c = (0, 0, 1)

