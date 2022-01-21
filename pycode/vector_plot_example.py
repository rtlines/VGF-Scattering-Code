# vector plot
import matplotlib.pyplot as plt
v = [1,1,0]



fig=plt.figure()
ax=plt.axes(projection = '3d')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
start = [0.0,0.0,0.0]
#ax.quiver(start[0],start[1],start[2],v[0],v[1],v[2])
ax.quiver(0,0,0,1,1,1)
ax.quiver(start[0],start[1],start[2],v[0],v[1],v[2])