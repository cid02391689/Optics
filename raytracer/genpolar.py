import numpy as np
def rtrings(rmax=5., nrings=5, multi=6):
    yield (0, 0) #origin
    for i in range(nrings): # 0~nrings-1
        r = rmax * (i + 1) / nrings
        N_points= (i+1)*multi
        thetas = np.linspace(0, 2*np.pi, N_points, endpoint=False)
        for theta in thetas:    
            yield r, theta
            # yield r * np.cos(theta), r * np.sin(theta), 0
    
