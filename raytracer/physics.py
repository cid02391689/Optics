import numpy as np
def refract(direc, normal, n_1, n_2):
    k1 = np.array(direc)/np.linalg.norm(direc)
    normal = np.array(normal)/np.linalg.norm(normal)
    cos_theta_1 = -np.dot(k1, normal)
    sin_theta_2 = n_1 / n_2 * np.sqrt(1 - cos_theta_1**2)
    if sin_theta_2 > 1:
        return None
    cos_theta_2 = np.sqrt(1 - sin_theta_2**2)
    k2 = n_1 / n_2 * k1 + (n_1 / n_2 * cos_theta_1 - cos_theta_2) * normal
    return k2 / np.linalg.norm(k2)
