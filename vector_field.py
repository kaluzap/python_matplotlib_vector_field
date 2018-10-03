import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

c_tau = 1.0
c_delta = 5.0

def vector_field(i_a, i_b, i_c):
    i_U = -(1.0/c_tau)*(i_a*i_a - i_c)*(i_a - i_b)
    i_V = (1.0/c_delta)*(i_a - i_b)
    return i_U, i_V


def vector_field_directions(i_a, i_b, i_c):
    
    i_U = -(1.0/c_tau)*(i_a*i_a - i_c)*(i_a - i_b)
    i_V = (1.0/c_delta)*(i_a - i_b)
    
    for i in range(np.size(i_U, 0)):
        for j in range(np.size(i_U, 1)):
            
            norma = np.sqrt(i_U[i][j]*i_U[i][j] + i_V[i][j]*i_V[i][j])
            
            if norma != 0:
                i_U[i][j] /= norma
                i_V[i][j] /= norma
            else:
                i_U[i][j] = 0.0
                i_V[i][j] = 0.0
    return i_U, i_V


def vector_field_intensity(i_a, i_b, i_c):
    
    i_U = -(1.0/c_tau)*(i_a*i_a - i_c)*(i_a - i_b)
    i_V = (1.0/c_delta)*(i_a - i_b)
    
    norma = np.sqrt(i_U*i_U + i_V*i_V)
            
    return norma


n = 1.4
X, Y = np.mgrid[-n:n:0.1, -n:n:0.1]
Xi, Yi = np.mgrid[-n:n:0.03, -n:n:0.03]

U, V = vector_field_directions(X, Y, 1.0)
Z = vector_field_intensity(Xi, Yi, 1.0)

#unstable line
x0 = [-1.0, 1.0]
y0 = [-1.0, 1.0]

#stable line1
x1 = [-n, -1.0]
y1 = [-n, -1.0]

#stable line2
x2 = [1.0, n]
y2 = [1.0, n]




# ========================= 
plt.clf()
fig = plt.figure(1)
# =========================

plt.subplot(2, 1, 1)

plt.plot(x0, y0, linestyle = '--', color = 'blue')
plt.plot(x1, y1, color = 'blue')
plt.plot(x2, y2, color = 'blue')

plt.plot(-1.0,-1.0, 'bo')
plt.text(-0.975, -1.125, 'A', color = 'blue')

plt.plot(1.0,1.0, 'bo')
plt.text(0.925, 1.05, 'B', color = 'blue')

plt.quiver(X, Y, U, V, linewidth=0.1)

plt.xlabel('$w$')
plt.ylabel("$w'$")

plt.xlim(-n, n)
#plt.xticks(())
plt.ylim(-n, n)
#plt.yticks(())

plt.text(-1.8, 1.3, '(a)', fontsize=15, color = 'black')


plt.subplot(2, 1, 2)


pc = plt.pcolormesh(Xi, Yi, Z, cmap = cm.gray)

# plt.colorbar(pc)

plt.plot(x0, y0, linestyle = '--', color = 'blue')
plt.plot(x1, y1, color = 'blue')
plt.plot(x2, y2, color = 'blue')

plt.plot(-1.0,-1.0, 'bo')
plt.text(-0.975, -1.125, 'A', color = 'blue')

plt.plot(1.0,1.0, 'bo')
plt.text(0.925, 1.05, 'B', color = 'blue')





plt.text(-1.8, 1.3, '(b)', fontsize=15, color = 'black')


plt.xlabel('$w$')
plt.ylabel("$w'$")

plt.xlim(-n, n)
#plt.xticks(())
plt.ylim(-n, n)
#plt.yticks(())

# ========================= 
cax = fig.add_axes([0.92, 0.11, 0.03, 0.35])
cb = fig.colorbar(pc, cax, orientation='vertical')
# =========================

cax.tick_params(labelsize=8)



fig.set_figwidth(5)
fig.set_figheight(10)


plt.savefig('field.png', bbox_inches='tight')


plt.show()


