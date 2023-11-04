import numpy as np
import matplotlib.pyplot as plt

def fuel_mass_calculator(maneuver, engine, g, wet_mass, dry_mass, delta_v):

    #Isp of our Liquid Hydrogen Fuel
    Isp = 455

    fuel_mass = wet_mass - dry_mass

    delta_v_array = np.linspace(0, delta_v, 1000)
    fuel_consumed = np.zeros(1000)  # Initialize an array to store fuel consumed

    for i, d in enumerate(delta_v_array):
        if fuel_mass <= 3000:
            Isp = 900

            #Above step ^ If the fuel mass equals or drops below 3000kg, 
            #we will switch to the fuel obtained from Starman

        fuel_consumed[i] = wet_mass * (1 - np.exp(-d / (Isp * g)))
        fuel_mass = wet_mass - fuel_consumed[i] - dry_mass
        new_mass = wet_mass - fuel_consumed[-1]
            
    print(f"The fuel remaining after completing the {maneuver} is {fuel_mass}kg")

    print(f"The fuel consumed by the {maneuver} maneuver is {fuel_consumed[-1]}kg")

    print(f"The new mass after the {maneuver} maneuver is {new_mass}kg")

    return fuel_consumed, new_mass

#specific gravity values
g_earth = 9.81
g_mars = 3.72076
g_sun = 274

wet_mass = 63800
#dry mass subject to change
dry_mass = 3000
print(f"The dry mass of the spacecraft is equal to {dry_mass}kg")

starman_fuel_mass = 3000

#Calculating fuel consumed and remaining mass for each maneuver
fuel_cons1, mass1 = fuel_mass_calculator("Parking to Hyperbolic", "falcon heavy", g_earth, wet_mass, dry_mass, 3598.4)

fuel_cons2, mass2 = fuel_mass_calculator("Starman Orbit Injection", "falcon heavy", g_sun, (mass1 + starman_fuel_mass), dry_mass, 701.17)

fuel_cons3, mass3 = fuel_mass_calculator("Starman to Mars", "falcon heavy", g_sun, mass2, dry_mass, 2030.5)

fuel_cons4, mass4 = fuel_mass_calculator("Parking Orbit Injection", "falcon heavy", g_mars, mass3, dry_mass, 111.75)

fuel_cons5, mass5 = fuel_mass_calculator("Parking Orbit Adjustment", "falcon heavy", g_mars, mass4, dry_mass, 1891.6)

fuel_cons6, mass6 = fuel_mass_calculator("Escape Mars' SOI", "falcon heavy", g_mars, mass5, dry_mass, 2142.1)

fuel_cons7, mass7 = fuel_mass_calculator("Earth Parking Orbit Injection", "falcon heavy", g_mars, mass6, dry_mass, 182.73)


fuel_list = [fuel_cons1[-1], fuel_cons2[-1], fuel_cons3[-1], fuel_cons4[-1], fuel_cons5[-1], fuel_cons6[-1], fuel_cons7[-1]]
delta_v_list = [3598.4, 701.17, 2030.5, 111.75, 1891.6, 2142.1, 182.73]

values = [1,2,3,4,5,6,7]
colors = []
for value in values:
    if value == 1:
        colors.append('blue')
    elif 2 <= value <= 3:
        colors.append('green')
    elif 4 <= value <= 6:
        colors.append("red")
    else:
        colors.append('blue')


plt.scatter(delta_v_list, fuel_list, color=colors)


plt.xlabel("Delta V (km/s)")
plt.ylabel("Fuel Consumed (kg)")
plt.title("Fuel Consumed vs Delta V")
plt.show()
plt.savefig("Fuel23.png")

