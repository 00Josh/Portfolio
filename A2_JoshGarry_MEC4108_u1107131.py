# A2 MEC4108 S1 2022
# Joshua Garry u1107131 Python Executable Assignment

#  To install requirements copy and past to command line:
#       pip install numpy matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Define some global functions
C2K = 273.15  # (K)

def loading_bar(x, total):
    n = 50
    bar = round((x + 1) * n / total)
    pct = str(round(100 * x / total, 3))
    print("\r", end='')
    print("Solving [" + '=' * bar + ' ' * (n - bar) + '] ' + pct + '%', end='')







print('Question 1:\n\n')
"""
Assumptions:
    That the heat conduction in the cylinder is two dimensional along the axis and radius of the cylinder.
    That the thermal properties of the cylinder are the same in all directions.
    The Fourier number is greater than 0.2 so that one-term solutions are appropriate.
"""
# plt.figure()
# plt.imshow(mpimg.imread('./Sketch_Q1.jpg'))

# Material : AISI 302
r = 0.100/2     # (m)
l = 0.150       # (m)
L = l/2         # m
Ti = 250        # (C)
T_inf = 35      # (C)
t = 8 * 60  # (seconds)
h = 205         # (W/[m^2 * K])
k = 15.18       # (W/[m*K]) at 300 K
a = 3.91*10**-6 # (m^2/s) Thermal diffusivity For AISI 302 at 300K from Fund of Heat and Mass Trf 7th Table A1 ~P984
rho = 8055      # kg/m^3
C_p = 480       # J/[kg*K]


# For the surface when tau > 0.2 the approximate dimensionless temperature formula is:
def theta_surf(x, L, tau, Bi):
    Bi_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    lambda1_range = [0.3111, 0.4328, 0.5218, 0.5932, 0.6533, 0.7051, 0.7506, 0.7910, 0.8274, 0.8603]
    A1_range = [1.0161, 1.0311, 1.0450, 1.0580, 1.0701, 1.0814, 1.0918, 1.1016, 1.1107, 1.1191]

    # if Bi > Bi_range[-1] or Bi < Bi_range[0]: print('ERROR: Bi No out of interp range')

    lambda1 = np.interp(Bi, Bi_range, lambda1_range)
    A1 = np.interp(Bi, Bi_range, A1_range)

    return A1 * np.e**(-lambda1**2 * tau) * np.cos(lambda1 * x/L)


# For a Cylinder the tau > 0.2 the approximate dimensionless temperature formula is:
def theta_cyl(r, r0, tau, Bi):
    Bi_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    lambda1_range = [1.22558, 1.5995, 1.7887, 1.9081, 1.9898, 2.0490, 2.0937, 2.1286, 2.1566, 2.1795]
    A1_range = [1.2071, 1.3383, 1.4191, 1.4698, 1.5029, 1.5253, 1.5411, 1.5526, 1.5611, 1.5677]
    J0_range = [1, 0.7665, 0.7196, 0.671, 0.6201, 0.559]
    n_range = [0, 1, 1.1, 1.2, 1.3, 1.4]

    # if Bi > Bi_range[-1] or Bi < Bi_range[0]: print('ERROR: Bi No out of interp range')
    lambda1 = np.interp(Bi, Bi_range, lambda1_range)
    A1 = np.interp(Bi, Bi_range, A1_range)

    # unnecessary for this question as r=0 => J0=1
    n = lambda1 * r/r0
    if n > n_range[-1] or n < n_range[0]: print('ERROR: n in J0(n) out of interp range')
    J0 = np.interp(n, n_range, J0_range)

    return A1 * np.e**(-lambda1**2 * tau) * J0


print('For the centre of the surface')
tau_suf = a * t / L**2     # N.D
# 1/Bi = k / (h * L)
Bi_suf = (h * L) / k
theta_suf_0_t = theta_surf(0, L, tau_suf, Bi_suf)
print('Tau  :', tau_suf)
print('Bi   :', Bi_suf)
print('1/Bi :', Bi_suf**-1)
print('theta:', theta_suf_0_t)

print('\nFor the centre of the cylinder')
tau_cyl = a * t / r**2     # N.D
Bi_cyl = (h * r) / k
theta_cyl_0_t = theta_cyl(0, r, tau_cyl, Bi_cyl)
print('Tau  :', tau_cyl)
print('Bi   :', Bi_cyl)
print('1/Bi :', Bi_cyl**-1)
print('theta:', theta_cyl_0_t)

''' At the intersection of the two solids:
        (T(0,0,t) - T_inf)/(T_i - T_inf) = theta_wall(0, t) * theta_cyl(0, t)
'''
T_0_0_t = theta_suf_0_t * theta_cyl_0_t * (Ti - T_inf) + T_inf
print('\nThe temperature of the geometric centre of the cylinder after 8 min (C)\n', T_0_0_t)

T_L_0_t = theta_surf(L, L, tau_suf, Bi_suf) * theta_cyl_0_t * (Ti - T_inf) + T_inf
print('\nThe temperature of the surface at centre of the cylinder after 8 min (C)\n', T_L_0_t)

# For the total heat loss in the system
V = np.pi * r**2 * l           # m^3
m = rho * V                    # kg
Q_max = m * C_p * (Ti - T_inf) # J
'''The total heat loss can be described by:
        (Q/Q_max)_total = (Q/Q_max)_wall + (Q/Q_max)_cyl * [1 - (Q/Q_max)_cyl]
                        = X + Y * [1 - Y]
        Q = (X + Y * [1 - Y]) * Q_max_total
        (Q/Q_max)_wall = 1 - theta(0, t)_wall * sin(lambda1_wall)/lambda1_wall
        (Q/Q_max)_cyl = 1 - theta(0, t)_wall * J1(lambda1_cyl)/lambda1_cyl
'''
Bi_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
lambda1_suf_range = [0.3111, 0.4328, 0.5218, 0.5932, 0.6533, 0.7051, 0.7506, 0.7910, 0.8274, 0.8603]
lambda1_cyl_range = [1.22558, 1.5995, 1.7887, 1.9081, 1.9898, 2.0490, 2.0937, 2.1286, 2.1566, 2.1795]
lambda1_suf = np.interp(Bi_suf, Bi_range, lambda1_suf_range)
lambda1_cyl = np.interp(Bi_cyl, Bi_range, lambda1_cyl_range)

J1_range = [0, 0.4400, 0.4709, 0.4983, 0.5220, 0.5419]
n_range = [0, 1, 1.1, 1.2, 1.3, 1.4]

t_arr = np.arange(2*60 , 20*60, 2)
Q = []
for t_x in t_arr:
    tau_suf = a * t_x / L ** 2  # N.D
    tau_cyl = a * t_x / r ** 2  # N.D
    X = 1 - theta_surf(0, L, tau_suf, Bi_suf) * np.sin(lambda1_suf)/lambda1_suf
    J1 = np.interp(lambda1_cyl, n_range, J1_range)
    Y = 1 - theta_cyl(0, r, tau_cyl, Bi_cyl) * J1/lambda1_cyl
    Q.append((X + Y * (1 - Y)) * Q_max)

# to convert from heat energy left to heat output in kJ
Q = [(Q_max - Q_x)/1000 for Q_x in Q]

plt.plot(t_arr/60, Q)
plt.title('Q1 Total heat loss')
plt.xlabel('Time (Minutes)')
plt.ylabel('Heat Energy (kJ)')
plt.grid()











print('\n\nQuestion 2:\n\n')
"""
For simplicity I will refer to the first case as:
    when the right surface is defined by a convection coefficient and air temperature of 'h' and 'T_infinity',
and the second case as:
    when the right side surface is defined by a consent heat flux of 'q'.

Assumptions:
    That the k value of the material is constant with temperature
    That the radiation heat transfer is negligible
    That the plate can be approximated as infinitely long
    That no heat generation occurs inside the pate
    That the material has consistent heat transfer characteristics in all directions
    That the plate can be assumed to be perfectly symmetrical
"""
# plt.figure()
# plt.imshow(mpimg.imread('./Sketch_Q2.jpg'))

T_bsurf = 225                       # (C)
T_inf = 22                          # (C)
h = 75                              # (W/[m^2 * K])
a = 3.91*10**-6                     # (m^2/s) Thermal diffusivity For AISI 302 at 300K from Fund of Heat and Mass Trf 7th Table A1 ~P984
k = 15.18                           # (W/[m*K]) at 300 K
q = 3500                            # (W/m^2) used in second case
x = 0.050                           # (m)
y = 0.030                           # (m)
t_check = [ 60, 5*60, 10*60, 30*60] # (s)

# T is a two-dimensional matrix of size 4x6 to include the lower surface
# (Just for fun, you can make this larger for a more accurate result)
# dims = [4, 6]
# dims = [7, 11]
dims = [14, 22]
T = np.zeros(dims)                  # (C)

'''
Based on Fourier number (Fo) or time stepping constant
    Fo = a * delta_t / delta_x^2
with stable solutions in two dimensions (Fund of Heat and Mass Trf 7th eq. 5.83 ~P332) for an  when
    Fo <= 1/4
    a * delta_t / delta_x^2 <= 1/4
    delta_t <= delta_x^2 / 4a
'''
delta_x = x/(dims[1] - 1)                # (m)
delta_y = y/(dims[0] - 1)                # (m)

max_delta_t = delta_x**2 / (4 * a)       # (s)
print("Maximum delta_t    : ", max_delta_t)
delta_t = np.round(0.6 * max_delta_t, 2) # (s)
print('Chosen delta_t     : ', delta_t)

Fo = a * delta_t / delta_x**2 # (N.D)
print('Fourier Number (Fo):', Fo)

Bi = h * delta_x / k # (N.D)  eq. 5.86 (Fund of Heat and Mass Trf 7th, ~P333)
print('Finite Biot No (Bi):', Bi)

print('\nStability criterion:')# will return a True/False
print('Fo <= 1/4          :', Fo <= 1/4)          # eq 5.83 Fund of Heat and Mass 7th Tbl 5.3 ~P334
print('Fo(3 + Bi) <= 3/4  :', Fo*(3 + Bi) <= 3/4) # eq 5.89
print('Fo(2 + Bi) <= 1/2  :', Fo*(2 + Bi) <= 1/2) # eq 5.91
print('Fo(1 + Bi) <= 1/4  :', Fo*(1 + Bi) <= 1/4) # eq 5.93


# Make a lookup list of important (m, n) matrix locations that will have unique T() functions
# make a list of the [m, j] matrix locations of inside cells
internal_nodes = []
for m in range(1, T.shape[0] - 1):
    for n in range(1, T.shape[1] - 1):
        internal_nodes.append((m, n))

top_surface_nodes = [(0, n) for n in range(1, T.shape[1] - 1)]
sym_surface_nodes = [(m, 0) for m in range(1, T.shape[0] - 1)]
right_surface_nodes = [(m, T.shape[1] - 1,) for m in range(1, T.shape[0] - 1)]
top_left_node = (0, 0)
top_right_node = (0, T.shape[1] - 1)
bottom_right_node = (T.shape[0] - 1, T.shape[1] - 1)

'''
Explicit method from first principals
    A*rho*C_p*delta_x*(T[m,n]_i+1 - T[m,n]_i)/delta_t = sum(Q_dot_i) + E_dot_gen

Because there is no heat generation in the material E_dot_gen = 0
    A*rho*V*C_p*(T[m,n]_i+1 - T[m,n]_i)/delta_t = sum(Q_dot_i)
    
Normalising both surfaces areas by depth
    rho * delta_x * C_p (T[m,n]_i+1 - T[m,n]_i)/delta_t = sum(q_dot_i)
    
    delta_x = delta_y
    a = k/(C_p * rho)
    Fo = a*t/delta_x^2                   (Fourier Number)
    Fo = k*t/(C_p * rho * delta_x^2)

    k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) = sum(q_dot_i)
    
    
For internal nodes:
    k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) = q_dot_top + q_dot_bottom + q_dot_left + q_dot_right
    k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) =  k*(T[m-1,n] - T[m,n])/delta_y
                                                + k*(T[m+1,n] - T[m,n])/delta_y
                                                + k*(T[m,n-1] - T[m,n])/delta_x
                                                + k*(T[m,n+1] - T[m,n])/delta_x
                                    
    T[m,n]_i+1 = Fo(T[m-1,n] + T[m+1, n] + T[m, n-1] + T[m, n+1] - 4*T[m,n]) + T[m, n]
    T[m,n]_i+1 = Fo(T[m-1,n] + T[m+1, n] + T[m, n-1] + T[m, n+1]) + (1 - 4*Fo)*T[m,n]
    
    
For top surface nodes:
    k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) =  h*(T_inf - T[m,n])
                                                + k*(T[m+1,n] - T[m,n])/delta_y
                                                + k*(T[m,n-1] - T[m,n])/(2 * delta_x)
                                                + k*(T[m,n+1] - T[m,n])/(2 * delta_x)
                                    
    Bi = delta_y * h / k (Biot No)
    
    T[m,n]_i+1 = Fo*( Bi*T_inf + T[m+1, n] + T[m, n-1]/2 + T[m, n+1]/2 - (Bi + 2)*T[m,n]) + T[m, n]
    T[m,n]_i+1 = Fo*( Bi*T_inf + T[m+1, n] + T[m, n-1]/2 + T[m, n+1]/2) + (1 - Bi*Fo - 2*Fo)*T[m,n]
    
    
For the right side surface nodes:
    Case 1
        The top surface node formula can be adjusted
        T[m,n]_i+1 = Fo*( T[m-1,n]/2 + T[m+1, n]/2 + T[m, n-1] + Bi*T_inf) + (1 - Bi*Fo -2*Fo)*T[m,n]
    
    Case 2
        k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) =  k*delta_x*(T[m-1,n] - T[m,n])/(2 * delta_y) 
                                                    + k*delta_x*(T[m+1,n] - T[m,n])/(2 * delta_y)  
                                                    + k*delta_y*(T[m,n-1] - T[m,n])/delta_x 
                                                    + q_dot_right
                
        T[m,n]_i+1 =   Fo * ( (T[m-1,n] - T[m,n])/2  
                            + (T[m+1,n] - T[m,n])/2 
                            + (T[m,n-1] - T[m,n])
                            q_dot_right / k ) + T[m,n]
                            
        T[m,n]_i+1 = Fo * (T[m-1,n]/2 + T[m+1,n]/2 + T[m,n-1] + q_dot_right / k) -2*Fo*T[m,n] + T[m,n]
        T[m,n]_i+1 = Fo * (T[m-1,n]/2 + T[m+1,n]/2 + T[m,n-1] + q_dot_right / k) + (1 -2*Fo) * T[m,n]


For the symmetrical surface cells:
    The formula for the internal nodes can be adjusted to 
    T[m,n]_i+1 = Fo(T[m-1,n] + T[m+1, n] + 2*T[m, n+1]) + (1 - 4*Fo)*T[m,n]


For the top left node:
    the top surface node formula can be adjusted to mirror the symmetry line
    T[m,n]_i+1 = Fo( Bi*T_inf + T[m+1, n] + T[m, n+1]) + (1 - Bi*Fo - 2*Fo)*T[m,n]


For the top right node:
    Case 1:
        k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) =  h*delta_y*(T_inf - T[m,n])/2
                                                    + k*delta_x*(T[m+1,n] - T[m,n])/(2 * delta_y)
                                                    + k*delta_y*(T[m,n-1] - T[m,n])/(2 * delta_x)
                                                    + h*delta_x*(T_inf - T[m,n])/2
                                        
        T[m,n]_i+1 = Fo * ( Bi*T_inf + T[m+1,n] + T[m,n-1] + Bi*T_inf - (2 -2*Bi)*T[m, n] )/2 + T[m,n]
        T[m,n]_i+1 = Fo * ( Bi*T_inf + T[m+1,n] + T[m,n-1] + Bi*T_inf)/2 + (1 - Fo + Bi*Fo)*T[m, n]
        
    Case 2:
        k * (T[m,n]_i+1 - T[m,n]_i)/(Fo * delta_x) =  h*delta_y*(T_inf - T[m,n])/2
                                                    + k*delta_x*(T[m+1,n] - T[m,n])/(2 * delta_y)
                                                    + k*delta_y*(T[m,n-1] - T[m,n])/(2 * delta_x)
                                                    + q*delta_x/2
                                        
        T[m,n]_i+1 = Fo * ( Bi*T_inf/2 + T[m+1,n]/2 + T[m,n-1]/2 + q*delta_x/(2*k) - 3*T[m, n] ) + T[m, n]
        T[m,n]_i+1 = Fo * ( Bi*T_inf/2 + T[m+1,n]/2 + T[m,n-1]/2 + q*delta_x/(2*k) ) + (1 - 3*Fo)*T[m, n]

'''

Q2_T_ANS = [[], []]

# Solve for scenario case 1 and 2 of this problem
for case in [0, 1]:
    print('\nCase:', case+1)
    ans_count = 0

    T[:, :] = T_bsurf  # set the temp of all the nodes in T to 255 C
    T2 = T.copy()

    # set the time range array
    t_range = np.arange(0, t_check[-1] + delta_t, delta_t)
    t_prev = 0
    for t in t_range:
        #  this will iterate though all (m, n) points and use the formula for the type of cell it is
        """
        Using The finite difference equations with time stepping:
        """
        for m in range(0, T.shape[0]):
            for n in range(0, T.shape[1]):
                if (m, n) in internal_nodes:
                    # Cells in the middle of the plate can use eq. 5.3 (Fund of Heat and Mass Trf 7th, ~P334)
                    # because delta_x = delta_y
                    T2[m, n] = Fo * (T[m-1, n] + T[m+1, n] + T[m, n-1] + T[m, n+1]) + (1 - 4*Fo) * T[m, n]

                elif (m, n) in sym_surface_nodes:
                    # Identical to the above formula with the   v   exception of the using the same point for both
                    # sides of the symmetry line
                    T2[m, n] = Fo * (T[m-1, n] + T[m+1, n] + 2 * T[m, n+1]) + (1 - 4*Fo) * T[m, n]

                elif (m, n) in right_surface_nodes:
                    if case == 0:
                        T2[m, n] = Fo*( T[m-1, n] + T[m+1, n] + 2*T[m, n-1] + 2*Bi*T_inf) + (1 - 2*Bi*Fo - 4*Fo) * T[m, n]

                    elif case == 1:
                        T2[m, n] = Fo * (T[m-1, n]/2 + T[m+1, n]/2 + T[m, n-1] + q/k) + (1 - 2*Fo) * T[m, n]

                elif (m, n) in top_surface_nodes:
                    T2[m, n] = Fo*( 2*Bi*T_inf + 2*T[m+1, n] + T[m, n-1] + T[m, n+1]) + (1 - 2*Bi*Fo - 4*Fo)*T[m, n]

                elif (m, n) == top_left_node:
                    # from eq 5.90 Fund of Heat and Mass Trf 7th ~P334
                    # because of symmetry T[m, n-1] = T[m, n+1], rearranging to:
                    T2[m, n] = Fo * ( 2*Bi*T_inf + 2*T[m+1, n] + 2*T[m, n+1]) + (1 - 2*Bi*Fo - 4*Fo)*T[m, n]

                elif (m, n) == top_right_node:
                    if case == 0:
                        T2[m, n] = 2*Fo * ( 2*Bi*T_inf + T[m+1, n] + T[m, n-1]) + (1 - 4*Fo - 4*Bi*Fo)*T[m, n]

                    elif case == 1:
                        T[m, n] = Fo * ( Bi*T_inf/2 + T[m+1, n]/2 + T[m, n-1]/2 + q/(2*k)) + (1 - 3*Fo)*T[m, n]

                loading_bar(t, t_range[-1])

        T = T2.copy()

        # save the T matrix at 1, 5, 10 and 30 min
        if t_prev < t_check[ans_count] <= t:
            ans_count += 1
            Q2_T_ANS[case].append(T.copy())

        t_prev = t


# Display the answer in the printout and plot the heat maps for case 1
for time, T_ans in zip(t_check, Q2_T_ANS[0]):
    print('\n\nFinal temperature matrix (C) at t =', time, 'for case 1')
    for i, temp in enumerate(T_ans.reshape([1, -1])[0, :18]):
        print('    T_' + str(i+1), ' = ', round(temp, 3))

    plt.matshow( T_ans, cmap='coolwarm')
    cbar = plt.colorbar()
    cbar.set_label('Temperature (C)') # , rotation=270
    plt.title('Temperature heatmap at ' + str(time) + ' seconds case 1')
    plt.xlabel('X Axis Node Index')
    plt.ylabel('Y Axis Node Index')

# Display the heat map for second case scenario at 10min
plt.matshow(Q2_T_ANS[1][2], cmap='coolwarm')
cbar = plt.colorbar()
cbar.set_label('Temperature (C)')  # , rotation=270
plt.title('Temperature heatmap at ' + str(10*60) + ' seconds case 2')
plt.xlabel('X Axis Node Index')
plt.ylabel('Y Axis Node Index')

# Print out the maximum node temp in case 2
T_c2 = Q2_T_ANS[1][2]
T_max_case2 = np.max(T_c2[:-1, :])
T_max_idxs = np.where(T_c2[:-1, :] == T_max_case2)
T_max_idx = T.shape[1] * T_max_idxs[0][0] + T_max_idxs[1][0]+1
print('\nCase 2 Maximum temperature cell:', round(T_max_case2, 2), 'C at T' + str(T_max_idx))









print('\n\nQuestion 3:\n\n')

"""
Assumptions:
    That there is no radiation present
    That the materials have constant properties
    That the system is in steady state
    That the ends of the plate can be idealised to be a 1
    That the plate is thin enough that the surface temperature can be assumed the same on both sides
"""
# plt.figure()
# plt.imshow(mpimg.imread('./Sketch_Q3.jpg'))

Q_dot = 2000 # (W)
L = 0.400    # (m)
h = 0.400    # (m)
t = 0.005    # (m)
T_inf = 25   # (C)
g = 9.81     # (m*s^-2)
A = l*h      # (m^2)

# from table A5 in Fund of heat and Mass Trasf assuming a surface temp of 350 (K)
k = 138*10**-3   # (W/[kg*K])
nu = 41.7*10**-6 # (W/[m*K])
Pr = 546         # N.D

# Iteratively solve for the vertical plate surface temp
T_s = 350 - C2K # (C)
for i in range(0, 10):
    T_f = (T_s + T_inf)/2                                                  # (C)
    beta = (T_s + C2K)**-1                                                 # (1/K)
    Gr = g * beta * (T_s - T_inf) * L**3 / nu**2                           # N.D (Grashof No)
    Ra = Gr*Pr                                                             # N.D (Rayleigh no)
    Nu = ( 0.825 + 0.387*Ra**(1/6)/((1 + (0.492/Pr)**(9/16))**(8/27) ))**2 # N.D (Nusselt no)
    # Nu = h*L_c/k
    h = Nu*k/L

    # Q/2 = A*h*(Ti - T_inf)
    T_s = Q_dot/(2*A*h) + T_inf

print('Vertical')
print('Grashof No                        :', Gr)
print('Rayleigh No                       :', Ra)
print('Nusselt No                        :', Nu)
print('Convection Coefficient (W/[m^2*K]):', h)
print('Surface Temp (C)                  :', T_s)


# Iteratively solve for the horizontal plate surface temp
T_s = 350 - C2K # (C)
for i in range(0, 10):
    T_f = (T_s + T_inf)/2                                                  # (C)
    beta = (T_s + C2K)**-1                                                 # (1/K)
    Gr = g * beta * (T_s - T_inf) * L**3 / nu**2                           # N.D (Grashof No)
    Ra = Gr*Pr                                                             # N.D (Rayleigh no)
    if Ra < 10*7:
        Nu_upper = 0.54*Ra**0.25                                           # N.D (Nusselt no)
    else:
        Nu_upper = 0.15*Ra**0.33333
    Nu_lower = 0.27 * Ra ** 0.25

    # Nu = h*L_c/k
    h_upper = Nu_upper * k / L
    h_lower = Nu_lower * k / L

    # Q = A*(h1 + h2)*(Ti - T_inf)
    T_s = Q_dot/(A*(h_upper + h_lower)) + T_inf

print('\nHorizontal')
print('Grashof No                        :', Gr)
print('Rayleigh No                       :', Ra)
print('Nusselt No')
print('    Top                           :', Nu_upper)
print('    Bottom                        :', Nu_lower)
print('Convection Coefficient (W/[m^2*K]):', h)
print('Surface Temp (C)                  :', T_s)

plt.show()
input('')
