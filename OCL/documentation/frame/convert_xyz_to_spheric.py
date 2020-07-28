import numpy as np
import math

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

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotated_vector(vector, axis, theta):
    return np.dot(rotation_matrix(axis,theta), vector)


# mm; Distance of global model "centre" to centre of ball/OSCAR
dx_from_center = 267.655

# load dataset
data = np.loadtxt("Koord_ball.csv")

# pick gobal xyz coordinates
xyz = data[:,4:]
# transform such that center of the ball is the origin
xyz[:,0] -= dx_from_center

# Rotate the axes such that it fits to our setup
v = np.array([[1,0,0],[2,0,1]])

xaxis  = [1, 0 , 0]
xtheta = np.pi/10. # rotation angle

yaxis = [0, 1, 0]
ytheta = -np.pi/2. # rotation angle

xyz_org = np.copy(xyz)
v_dummy = np.zeros((1,3))
# print v_dummy
for i, vec in enumerate(xyz_org):
    v_dummy = rotated_vector(vec,xaxis,xtheta)
    xyz[i] = rotated_vector(v_dummy,yaxis,ytheta)

sph_coords = append_spherical_np(xyz)
# sph_coords[:,4:] /= np.pi

# np.set_printoptions(precision=2,suppress=True)
np.set_printoptions(suppress=True)

#convert to degree
sph_coords[:,-2:] *= 180/np.pi
# print sph_coords[:,-2:] # print (theta,phi)
# print sph_coords

n_hexa = 20
n_pent = 12
n_tot  = n_pent + n_hexa
for i in range(n_hexa):
	print("frameHexagon_theta[{:2d}]\t= {:.5f}*deg;\
			\nframeHexagon_phi[{:2d}]\t= {:.5f}*deg;\n\
			 ".format(i,sph_coords[i,-2],i,sph_coords[i,-1],i,0))

for n in range(n_pent):
	i = n_tot - n_pent + n
	print("framePentagon_theta[{:2d}]\t= {:.5f}*deg;\
			\nframePentagon_phi[{:2d}]\t= {:.5f}*deg;\n\
			 ".format(n,sph_coords[i,-2],n,sph_coords[i,-1],n,0))

for n in range(n_tot):
    print("OCLLaBr3_presence[{:2d}]\t\t= true;\
             \nOCLCollimator_presence[{:2d}]\t= true;\
             \nOCLLaBr3_Distance[{:2d}]\t\t= {};\
             \nOCLLaBr3_theta[{:2d}]\t\t\t= {:.5f}*deg;\
             \nOCLLaBr3_phi[{:2d}]\t\t\t= {:.5f}*deg;\n\
             ".format(n,n,n,"20*cm",n,sph_coords[n,-2],n,sph_coords[n,-1]))
