# This code was developed for the GNSS Coursework in the Space Data Systems and Processing 
# It comprises a basic iterative algorithm to calculate the position of a GNSS user based on satellite and ground data

# imports
import numpy as np

# define initial data structures
pos_est = np.zeros(4) 
tolerance = 1e-10  # tolerance for convergence
max_iter = 1000  # max number of iterations
iter = 0  

# initial data
sat_pos = np.array([[10053355.1020386, 15056487.9624351, 25504822.96],
                    [20655550.9892971, 5147316.41599376, 16608524.31],
                    [21764223.9221401, -6594586.14401236, 15678237.26],
                    [10509172.2642949, -14025870.89, 25122346.95]])

range_meas = np.array([26226554.8, 21042234.5912586, 21770925.1002535, 25287580.1140081])



while iter < max_iter:
    # pseudo-range estimates
    range_est = np.zeros(4)
    for s in range(len(range_est)):
        range_est[s] = np.sqrt((sat_pos[s][0] - pos_est[0])**2 + 
                               (sat_pos[s][1] - pos_est[1])**2 + 
                               (sat_pos[s][2] - pos_est[2])**2) + pos_est[3]

    # range difference
    delta_range = np.subtract(range_meas, range_est)

    # measurement matrix
    H = np.zeros((4, 4))
    for s in range(4):
        for i in range(3):
            H[s][i] = -(sat_pos[s][i] - pos_est[i]) / range_est[s]
        H[s][3] = 1  # last column for clock bias

    # inverse measurement matrix
    H_inv = np.linalg.inv(H)

    # update position estimate
    new_pos_est = pos_est.copy()  
    new_pos_est += np.dot(H_inv, delta_range)

    # check for convergence
    if np.linalg.norm(new_pos_est - pos_est) < tolerance:
        break

    pos_est = new_pos_est 
    iter += 1 


# convert to curvilinear coordinates
print(pos_est)
a = 6378137
e = 0.081819191
c = np.sqrt(1 - e**2)

zeta = np.arctan(pos_est[2] / (c * np.sqrt(pos_est[0]**2 + pos_est[1]**2)))
easin = e**2 * a * (np.sin(zeta))**3
eacos = e**2 * a * (np.cos(zeta))**3

lat = np.arctan((c * pos_est[2] + easin) / (c * (np.sqrt(pos_est[0]**2 + pos_est[1]**2) - eacos)))
long = np.arctan2(pos_est[1], pos_est[0])

R = a / np.sqrt(1 - e**2 * (np.sin(lat)**2))

height = np.sqrt(pos_est[0]**2 + pos_est[1]**2) / np.cos(lat) - R

pos_est_curv = np.array([(180 / np.pi) * lat, (180 / np.pi) * long, height])



# Now print the results
print("Final position estimate:", pos_est)
print("Number of iterations:", iter)
print("Final delta_range:", delta_range)
print("Final measurement matrix H:\n", H)
print("Final inverse measurement matrix H_inv:\n", H_inv)
print("Final position estimate (curvilinear):", pos_est_curv)