import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

a = 2
lx = 2.0
ly = 1.0

try:
    T = float(input("Enter the total simulation time T: "))
    if T <= 0:
        raise ValueError("T must be positive.")

    Nx = int(input("Enter the number of grid points in x-direction: "))
    Ny = int(input("Enter the number of grid points in y-direction: "))
    Nt = int(input("Enter the number of time steps Nt: "))
    if Nx < 3 or Ny < 3 or Nt < 1:
        raise ValueError("Nx and Ny must be at least 3, Nt at least 1.")

except Exception as e:
    print("❌ Invalid input:", e)
    exit()

hx = lx / (Nx - 1)
hy = ly / (Ny - 1)
ht = T / Nt

print("hx = ", hx)
print("hy = ", hy)
print("ht = ", ht)

# Stability check
r = a * ht * np.sqrt(1/hx**2 + 1/hy**2)
print(f"Stability factor r = {r:.4f}")
if r > 1:
    raise SystemExit("❌ Warning: The scheme is unstable! Exiting...")

x = np.linspace(0, lx, Nx)
y = np.linspace(0, ly, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')

def g(x, y):
    return np.sin(np.pi*x/2) * np.sin(np.pi*y/2)
def q(x, y):
    return x * y

def d2_g(x, y):
  return -((np.pi / 2) ** 2) * np.sin(np.pi * x / 2) * np.sin(np.pi * y / 2)

def myu1(y, t):
  return 0
def myu2(y, t):
  return 2 * t * y
def myu3(x, t):
  return 0
def myu4(x, t):
  return t * x

def exact_solution(x, y, t):
    return t*x*y + np.cos(np.pi*np.sqrt(2)*t) * np.sin(np.pi*x/2) * np.sin(np.pi*y/2)

#U_all = np.zeros((Nt + 1, Nx, Ny))
U_prev = np.zeros((Nx, Ny))
U = np.zeros((Nx, Ny))
U_next = np.zeros((Nx, Ny))

U_exact = np.zeros((Nx, Ny))

# Initial condition
for m in range(Nx):
    for n in range(Ny):
        U_prev[m, n] = g(x[m], y[n])
        # First order initial condition
        U[m, n] = U_prev[m, n] + ht * q(x[m], y[n])
        # Second order initial condition
        #U[m, n] = (U_prev[m, n] + ht * q(x[m], y[n]) + ((a * ht) ** 2 ) / 2 * (d2_g(x[m], y[n]) + d2_g(x[m], y[n])))

# U_all[0] = U_prev
# U_all[1] = U

for k in range(1, Nt + 1):
        t = k * ht
        U_next.fill(0)

        U_next[1:-1, 1:-1] = (2 * U[1:-1, 1:-1] - U_prev[1:-1, 1:-1] +
                          (a * ht)**2 * (
                              (U[2:, 1:-1] - 2 * U[1:-1, 1:-1] + U[:-2, 1:-1]) / hx**2 +
                              (U[1:-1, 2:] - 2 * U[1:-1, 1:-1] + U[1:-1, :-2]) / hy**2))

        if k != Nt:
           U_prev[1:-1, 1:-1] = U[1:-1, 1:-1]
           U[1:-1, 1:-1] = U_next[1:-1, 1:-1]
        else:
            U[1:-1, 1:-1] = U_next[1:-1, 1:-1]

        # Boundary conditions
        U[0, :] = myu1(y, t)
        U[-1, :] = myu2(y, t)
        U[:, 0] = myu3(x, t)
		# First order boundary condition
        U[:, -1] = U[:, -2] + hy * myu4(x, t)
        # Second order boundary condition
        #U[:, -1] = (4 * U[:, -2] - U[:, -3] + 2 * hy * myu4(x, t)) / 3

        # U_all[k] = U

U_exact = np.zeros((Nx, Ny))
for m in range(Nx):
    for n in range(Ny):
        U_exact[m, n] = exact_solution(x[m], y[n], t)

# Max error
max_error = np.max(U - U_exact)
min_error = np.min(U - U_exact)
# Print maximum error
print(f"\n✅ Max error at t = {t:.3f}: {max_error}")
print(f"\n✅ Min error at t = {t:.3f}: {min_error}")

fig = plt.figure(figsize=(18, 6))

ax1 = fig.add_subplot(131, projection='3d')
ax1.plot_surface(X, Y, U, cmap='plasma')
ax1.set_title(f'Numerical Solution at t = {T}')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('u')

ax2 = fig.add_subplot(132, projection='3d')
ax2.plot_surface(X, Y, U_exact, cmap='viridis')
ax2.set_title(f'Exact Solution at t = {T}')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('u_exact')

ax3 = fig.add_subplot(133, projection='3d')
error = (U - U_exact)
ax3.plot_surface(X, Y, error, cmap='hot')
ax3.set_title(f'Error at t = {T}')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_zlabel('Error')

plt.tight_layout()
plt.show()