import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
from scipy.optimize import minimize
import FindStateVector as fsv
import PorkchopPlot as pp

# Define combined objective function
def total_delta_v(departure_dates):
    # Extract individual departure dates for Earth to Starman, Starman to Mars, Mars to Earth
    es_date = departure_dates[0]
    sm_date = departure_dates[1]
    me_date = sm_date + dt.timedelta(days = 90)  # 3 months later

    # Calculate delta-v for each section
    es_dv = calc_dv('Earth', 'Starman', es_date)
    sm_dv = calc_dv('Starman', 'Mars', sm_date)
    me_dv = calc_dv('Mars', 'Earth', me_date)

    # Return sum of delta-v for all sections
    return es_dv + sm_dv + me_dv

# Helper function to calculate delta-v for a given section
def calc_dv(origin, destination, departure_date):
    # Set up PorkchopPlot for the specific section
    sp_porkchop = pp.PorkchopPlot(origin, destination, departure_date.year, departure_date.month, departure_date.day,departure_date.year, departure_date.month, departure_date.day + 1, departure_date.year + 1, 120, 15, 24 * 30, 60, "pro", 120)

    # Get the minimum delta-v for the specific section
    sp_porkchop.get_min_dv()
    min_dv = sp_porkchop.min_dv

    return min_dv

# 3-4 month stay time constraint
def constraint_mars_stay(departure_dates):
    starman_to_mars_date = departure_dates[1]
    mars_to_earth_date = starman_to_mars_date + dt.timedelta(days = 90)  # 3 months later

    return (mars_to_earth_date - starman_to_mars_date).days - 90  # Constraint: 3-4 months separation

# Initial guess for departure dates
initial_departure_dates = [dt.date(2029, 1, 1), dt.date(2031, 1, 1), dt.date(2033, 1, 1)]

# Set bounds
bounds = [(dt.date(2029, 1, 1), dt.date(2033, 1, 1)),
          (dt.date(2030, 1, 1), dt.date(2035, 1, 1)),
          (dt.date(2032, 1, 1), dt.date(2040, 1, 1))]

# Define optimisation problem
optimization_result = minimize(total_delta_v, initial_departure_dates, method='SLSQP', bounds=bounds, constraints={'type': 'eq', 'fun': constraint_mars_stay})

# Extract optimised departure dates
optimized_departure_dates = optimization_result.x

print("Optimized Departure Dates:")
print("Earth to Starman:", optimized_departure_dates[0])
print("Starman to Mars:", optimized_departure_dates[1])
print("Mars to Earth:", optimized_departure_dates[1] + dt.timedelta(days = 90))