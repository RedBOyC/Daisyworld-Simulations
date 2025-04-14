import numpy as np
import math
import pandas as pd
from scipy.constants import sigma
import matplotlib.pyplot as plt

"""##constants"""

#Constants

sigma=sigma  #kg s^-3 K^-4
S=917  #W/m^2, solar constant
Aw=0.75  #albedo of white daisies
Ag=0.5   #albedo of bare ground
Ab=0.25  #albedo of black daisies
q_prime=20 # helps to calculate the local temperature as funct of local albedo
P=1.0 # proportion of the area with fertile ground
gamma=0.3 # death rate per unit time
T_opt = 22.5 # Optimal growth temperature [°C]
T_min = 5.0 # Minimum viable temperature [°C]
T_max = 40.0 # Maximum viable temperature [°C]

"""##Eqns"""

#Define a bunch of functions to make the code easier to read later

def growth_rate(T_L):  #eqn 3
  if T_min<= T_L <=T_max:
    return 1 - 0.003265*(22.5-T_L)**2
  else:
    return 0    ##if eqn 3 is less than 0, we can't have growth be <0

def albedo(aw, ab):  #eqn 5
  ag=1-aw-ab
  planetary_albedo=ag*Ag+ aw*Aw+ ab*Ab
  return planetary_albedo

def area_fertile(aw, ab):  #eqn 2
  area1=P -aw- ab
  return area1

#eqn 1:
def rate_black_daisies(ab, x, B):  #get x using eqn2 B from eqn 3
  rate_of_black_daisies=ab*((x*B) - gamma)
  return rate_of_black_daisies

def rate_white_daisies(aw, x, B): #get x using eqn2 B from eqn 3
  rate_of_white_daisies=aw*((x*B) - gamma)
  return rate_of_white_daisies

def effective_temp(A, L):  #eqn 4, get A from eqn 5
  temp_effective= ((S*L*(1-A))/sigma)**0.25 -273 #[°C]
  return temp_effective

def local_temp( A, A_i, T_e): #eqn 7, A_i is Ab or Aw
  T_i=q_prime*(A-A_i)+T_e #T_i is local temp of daisy i where i is b or w
  return T_i

"""## Simulation (guiding qn3)"""

#Only black daisies are present
L_values = np.arange(0.5, 2.05, 0.05)
final_temp = []
final_area = []

for L in L_values:
    # 1% black, 0% white
    ab = 0.01
    aw = 0

    # Track previous value to detect convergence
    prev_ab = 0

    for t in range(1000):
        A = albedo(aw, ab) #planet albedo
        T_e = effective_temp(A, L)  #effective temp at this L
        T_b = local_temp(A, Ab, T_e) #local temp of black daisies
        B = growth_rate(T_b)  #growth rate of black daisy
        x = area_fertile(aw, ab)   #fertile area available
        dab=rate_black_daisies(ab, x, B)*0.9  #dt=0.9
        ab+=dab #Euler method

        # Convergence test
        if abs(ab - prev_ab) == 0:
            break
        prev_ab = ab

    #Final steady state for this L
    A_final = albedo(aw, ab)
    T_e_final = effective_temp(A_final, L)

    #Record results
    final_area.append(ab * 100)  # Convert to %
    final_temp.append(T_e_final)


# Plot results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

# area coverage
ax1.plot(L_values, final_area, 'k-', label='BLACK')
ax1.set_ylabel('AREA (%)')
ax1.set_ylim(0, 70)
ax1.grid(True)
ax1.legend()

# temperature
ax2.plot(L_values, final_temp, 'k-')
ax2.set_xlabel('SOLAR LUMINOSITY')
ax2.set_ylabel('TEMP. [°C]')
ax2.set_ylim(0, 60)
ax2.grid(True)

plt.tight_layout()
plt.show()

#Only white daisy are available
# Increasing luminosity
L_values_inc = np.arange(0.5, 2.05, 0.05)
temp_inc = []
area_inc = []

# Decreasing luminosity
L_values_dec = np.arange(2.00, 0.45, -0.05)
temp_dec = []
area_dec = []

# First part: run with increasing luminosity, 1% white, 0% black
for L in L_values_inc:
  # 1% white, 0% black
  aw = 0.01
  ab = 0

  # Track previous value to detect convergence
  prev_aw = 0

  for t in range(1000):
    A = albedo(aw, ab)  # Planet albedo
    T_e = effective_temp(A, L)  # Effective temp at this L
    T_w = local_temp(A, Aw, T_e)  # Local temp of white daisies
    B = growth_rate(T_w)  #growth rate for white daisy
    x = area_fertile(aw, ab)  # Area available
    daw = rate_white_daisies(aw, x, B) * 0.9  # dt=0.9
    aw += daw  # Euler method

    # Convergence test
    if abs(aw - prev_aw) == 0:
      break
    prev_aw = aw


  # Record results
  A_final = albedo(aw, ab)
  T_e_final = effective_temp(A_final, L)

  area_inc.append(aw * 100)
  temp_inc.append(T_e_final)

# Now run with decreasing luminosity (2nd part)
# Start with the final white daisy population from increasing run
last_aw = aw

for L in L_values_dec:
  # Use final population from increasing run as starting point
  aw = last_aw
  ab = 0

  # Track previous value to detect convergence
  prev_aw = 0

  for t in range(1000):
    A = albedo(aw, ab)  # Planet albedo
    T_e = effective_temp(A, L)  # Effective temp at this L
    T_w = local_temp(A, Aw, T_e)  # Local temp of white daisies
    B = growth_rate(T_w)  #growth rate for white daisy
    x = area_fertile(aw, ab)  # fertile area available
    daw = rate_white_daisies(aw, x, B) * 0.9  # dt=0.9
    aw += daw  # Euler method

    # Convergence test
    if abs(aw - prev_aw) == 0:
      break
    prev_aw = aw


  #Final steady state for this L
  A_final = albedo(aw, ab)
  T_e_final = effective_temp(A_final, L)

  # Record results
  area_dec.append(aw * 100)
  temp_dec.append(T_e_final)

  # Update last_aw for the next iteration
  last_aw = aw

# Plot results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

# Plot area coverage
ax1.plot(L_values_inc, area_inc, 'k-', label='WHITE (increasing)')
ax1.plot(L_values_dec, area_dec, 'k--', label='WHITE (decreasing)')
ax1.set_ylabel('AREA (%)')
ax1.set_ylim(0, 70)
ax1.grid(True)
ax1.legend()

# Plot temperature
ax2.plot(L_values_inc, temp_inc, 'k-', label='Increasing')
ax2.plot(L_values_dec, temp_dec, 'k--', label='Decreasing')
ax2.set_xlabel('SOLAR LUMINOSITY')
ax2.set_ylabel('TEMP. [°C]')
ax2.set_ylim(0, 60)
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()

#Both black and white daisies are available
#Range of luminosity values
L_values=np.arange(0.5, 2.05, 0.05)

# Arrays to store results
temps = []
black_areas = []
white_areas = []

for L in L_values:

  # Initial conditions: 1% black, 1% white
  ab = 0.01
  aw = 0.01

  # Track previous values to detect convergence
  prev_ab = 0
  prev_aw=0

  for t in range(1000):
    A=albedo(aw, ab) #Planet albedo
    T_e=effective_temp(A, L) #Effective temp at this L
    T_b=local_temp(A, Ab, T_e) #Local temp of black daisies
    T_w=local_temp(A, Aw, T_e) #Local temp of white daisies
    B_b=growth_rate(T_b) #growth rate for black daisy
    B_w=growth_rate(T_w) #growth rate for white daisy
    x=area_fertile(aw, ab) #Fertile area available
    dab=rate_black_daisies(ab, x, B_b)*0.9 #dt=0.9
    daw=rate_white_daisies(aw, x, B_w)*0.9 #dt=0.9
    ab+=dab #Euler method
    aw+=daw #Euler method

    #Convergence test
    if abs(ab - prev_ab) == 0 and abs(aw - prev_aw) == 0:
      break
    prev_ab = ab
    prev_aw = aw

  #Final steady state for this L
  A_final=albedo(aw, ab)
  T_e_final=effective_temp(A_final, L)

  #Store results
  temps.append(T_e_final)
  black_areas.append(ab * 100) #Convert to percentage
  white_areas.append(aw * 100) #Convert to percentage

#Plot results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot daisy areas (top subplot)
ax1.plot(L_values, black_areas, 'k-', label='BLACK')
ax1.plot(L_values, white_areas, 'm-', label='WHITE')
ax1.set_ylabel('AREA (%)')
ax1.set_ylim(0, 80)
ax1.legend()
ax1.grid(True)

# Plot temperatures (bottom subplot)
ax2.plot(L_values, temps, 'k-')
ax2.set_xlabel('SOLAR LUMINOSITY')
ax2.set_ylabel('TEMP [°C]')
ax2.set_ylim(0, 70)
ax2.grid(True)

plt.tight_layout()
plt.show()
