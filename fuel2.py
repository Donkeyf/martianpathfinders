import numpy as np
import matplotlib.pyplot as plt

def fuel_mass_calculator(maneuver, engine, g, wet_mass, dry_mass, delta_v):

    Isp = 455

    fuel_mass = wet_mass - dry_mass

    delta_v_array = np.linspace(0, delta_v, 1000)
    fuel_consumed = np.zeros(1000)  # Initialize an array to store fuel consumed

    for i, d in enumerate(delta_v_array):
        if fuel_mass < 3000:
            Isp = 900

        fuel_consumed[i] = wet_mass * (1 - np.exp(-d / (Isp * g)))
        fuel_mass = wet_mass - fuel_consumed[i] - dry_mass
        new_mass = wet_mass - fuel_consumed[-1]
            
    print(fuel_mass)

    print(f"The fuel consumed by the {maneuver} maneuver is {fuel_consumed[-1]}kg")

    print(f"The new mass after the {maneuver} maneuver is {new_mass}kg")

    return fuel_consumed, new_mass

g_earth = 9.81
g_mars = 3.72076


#Spacecraft was launched from Earth, hence we use g_earth
fuel_cons1, mass1 = fuel_mass_calculator("Earth to Starman", "falcon heavy", g_earth, 63800, 1200, 3970.553988)

fuel_cons2, mass2 = fuel_mass_calculator("Starman to Mars", "falcon heavy", g_earth, (mass1 + 3000), 1200, 5429.575)

#Spacecraft launched from Mars orbit, hence we use g_mars
fuel_cons3, mass3 = fuel_mass_calculator("Mars to Earth", "falcon heavy", g_mars, mass2, 1200, 5935.492)

fuel_cons_list = [fuel_cons1[-1], fuel_cons2[-1], fuel_cons3[-1]]
delta_v_list = [3970.553988, 5429.575, 5935.492]
plt.figure()
plt.scatter(delta_v_list, fuel_cons_list)
plt.xlabel("Delta V (km/s)")
plt.ylabel("Fuel Consumed (kg)")


plt.savefig("fuel.png")


