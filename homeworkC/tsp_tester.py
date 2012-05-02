from math import *;

p_test = [[1,2],[2,1],[2,2]]



def distance(p):
    d = [[0 for j in range(len(p))] for i in range(len(p))]
    for i in range(len(p)):
        for j in range(len(p)):
            if(i!=j):
                d[i][j] = sqrt( (p[i][0]-p[j][0])**2 + (p[i][1]-p[j][1])**2)
            else:
                d[i][j] = 13371337;
    return d
    
def argmin(sequence, fn=None):
    """Two usage patterns:
    argmin([s0, s1, ...], fn)
    argmin([(fn(s0), s0), (fn(s1, s1), ...]) 
    Both return the si with lowest fn(si)"""
    if fn is None:
        return min(sequence)[1]
    else:
        return min((fn(e), e) for e in sequence)[1]

def nn(p,start): #only uses indces internally
    tot_dist = 0
    walk = []
    dist = distance(p)

    last_point = start
    walk.append(start)
    
    left_points = range(len(p))
    left_points.remove(start)
    for i in range(len(p)-1): 
        euclid_dist = lambda x: dist[last_point][x] 
        closest_point = argmin( left_points , euclid_dist ) #return index of the closest point
        tot_dist+=dist[last_point][closest_point] 
        last_point = closest_point
        walk.append(closest_point)
        left_points.remove(closest_point)

    tot_dist+=dist[last_point][start]
    walk.append(start)
    
    return (tot_dist,walk)

def opt(p):
    return p


print nn(p_test,0)

