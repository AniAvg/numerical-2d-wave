import numpy as np
import matplotlib.pyplot as plt

a = 2
lx = 2.0
ly = 1.0

T = float(input("Enter the total simulation time T: "))
Nx = int(input("Enter the number of grid points in x-direction: "))
Ny = int(input("Enter the number of grid points in y-direction: "))
Nt = int(input("Enter the number of time steps Nt: "))

hx = lx / (Nx - 1)
hy = ly / (Ny - 1)
ht = T / Nt

print(f"hx = {hx}")
print(f"hy = {hy}")
print(f"ht = {ht}")

# Stability check
r = a * ht * np.sqrt(1/hx**2 + 1/hy**2)
print(f"Stability factor r = {r:.4f}")
if r > 1:
   raise SystemExit("❌ Scheme is unstable!")

x = np.linspace(0, lx, Nx)
y = np.linspace(0, ly, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Initial and boundary condition functions
def g(x, y):
  return np.sin(np.pi * x / 2) * np.sin(np.pi * y / 2)
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
   return t * x * y + np.cos(np.pi * np.sqrt(2) * t) * np.sin(np.pi * x / 2) * np.sin(np.pi * y / 2)

U_all = np.zeros((Nt + 1, Nx, Ny))

# Initial condition
for i in range(Nx):
    for j in range(Ny):
        U_all[0, i, j] = g(x[i], y[j])

for i in range(Nx):
    for j in range(Ny):
        # First order initial condition
        U_all[1, i, j] = U_all[0, i, j] + ht * q(x[i], y[j])
        # Second order initial condition
        #U_all[1, i, j] = (U_all[0, i, j] + ht * q(x[i], y[j]) + ((a * ht) ** 2 ) / 2 * (d2_g(x[i], y[j]) + d2_g(x[i], y[j])))

for k in range(1, Nt):
    t = (k + 1) * ht

    U_all[k + 1, 1:-1, 1:-1] = (
        2 * U_all[k, 1:-1, 1:-1] - U_all[k - 1, 1:-1, 1:-1] +
        (a * ht) ** 2 * (
            (U_all[k, 2:, 1:-1] - 2 * U_all[k, 1:-1, 1:-1] + U_all[k, :-2, 1:-1]) / hx ** 2 +
            (U_all[k, 1:-1, 2:] - 2 * U_all[k, 1:-1, 1:-1] + U_all[k, 1:-1, :-2]) / hy ** 2
        )
    )

    # Boundary conditions
    U_all[k + 1, 0, :] = myu1(y, t)
    U_all[k + 1, -1, :] = myu2(y, t)
    U_all[k + 1, :, 0] = myu3(x, t)
    # First order boundary condition
    U_all[k + 1, :, -1] = U_all[k + 1, :, -2] + hy * myu4(x, t)
    # Second order boundary condition
    #U_all[k + 1, :, -1] = (4 * U_all[k + 1, :, -2] - U_all[k + 1, :, -3] + 2 * hy * myu4(x, t)) / 3

# Exact solution
U_exact = np.zeros((Nx, Ny))
t_final = Nt * ht
for i in range(Nx):
    for j in range(Ny):
        U_exact[i, j] = exact_solution(x[i], y[j], t_final)

# Max error
error = np.abs(U_all[Nt] - U_exact)
print(f"\n✅ Max absolute error at t = {t_final:.3f}: {np.max(error)}")

fig = plt.figure(figsize=(18, 6))
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot_surface(X, Y, U_all[Nt], cmap='plasma')
ax1.set_title(f'Numerical Solution at t = {T}')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('u')

ax2 = fig.add_subplot(132, projection='3d')
ax2.plot_surface(X, Y, U_exact, cmap='viridis')
ax2.set_title('Exact Solution')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('u_exact')

ax3 = fig.add_subplot(133, projection='3d')
ax3.plot_surface(X, Y, error, cmap='hot')
ax3.set_title('Error')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_zlabel('Error')

plt.tight_layout()
plt.show()