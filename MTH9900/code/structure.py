# Data structure: path and index
def total_num_nodes(n):
    return (4**(n+1)-1)//3

def start_index_at_step(n):
    return total_num_nodes(n-1)

def index_from_path(path):
    n = len(path)
    base = (4**n-1)//3
    sub = 0
    for i in range(n):
        assert path[i] >=0 and path[i] < 4
        sub = sub*4 + path[i]
    return base+sub

def path_from_index(index):
    n = 0
    base = 0
    while True:
        n += 1
        next_base = (4**n-1)//3
        if next_base > index: break
        base = next_base
    n = n-1    
    sub = index - base
    path = []
    while len(path) < n:
        i = sub % 4
        path.append(i)
        sub = sub // 4
    return path[::-1]
