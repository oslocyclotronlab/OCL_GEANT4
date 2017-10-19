import numpy as np

def append_spherical_np(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    # ptsnew = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,4] = np.arctan2(np.sqrt(xy),xyz[:,2])  # theta (0,pi); from z axis down
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])  # phi (-pi,pi)

    # map phi (-pi,pi) to (0,2*pi) 
    gt_idx = ptsnew[:,-1] < 0	# get rows where collum [:,i] is negative
    ptsnew[gt_idx,-1] = 2.*np.pi + ptsnew[gt_idx,-1]

    return ptsnew

# mm; Distance of global model "centre" to centre of ball/OSCAR
dx_from_center = 267.655

# load dataset
data = np.loadtxt("Koord_ball.csv")

# pick gobal xyz coordinates 
xyz = data[:,4:]
# transform such that center of the ball is the origin
xyz[:,0] -= dx_from_center

sph_coords = append_spherical_np(xyz)
# sph_coords[:,4:] /= np.pi

# np.set_printoptions(precision=2,suppress=True)
np.set_printoptions(suppress=True)

#convert to degree 
sph_coords[:,-2:] *= 180/np.pi
# print sph_coords[:,-2:] # print (theta,phi)
print sph_coords



# n_hexa = 20
# n_pent = 12
# for i in range(n_hexa):
# 	print    "frameHexagon_theta[%2i]\t= %f*deg;\
# 			\nframeHexagon_phi[%2i]\t= %f*deg;\
# 			\nframeHexagon_psi[%2i]\t= %i/twopi;\n\
# 			 " % (i,sph_coords[i,-2],i,sph_coords[i,-1],i,0)

# for n in range(n_pent):
# 	i = len(sph_coords) - n_pent + n
# 	print    "framePentagon_theta[%2i]\t= %f*deg;\
# 			\nframePentagon_phi[%2i]\t= %f*deg;\
# 			\nframePentagon_psi[%2i]\t= %i/twopi;\n\
# 			 " % (n,sph_coords[i,-2],n,sph_coords[i,-1],n,0)