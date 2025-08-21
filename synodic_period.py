def synodic_period_inferior(T_planet, T_earth):
    """
    Calculate synodic period for an inferior planet.
    T_planet: orbital period of planet (days)
    T_earth: orbital period of Earth (days)
    Returns: synodic period (days)
    """
    return 1 / (1 / T_planet - 1 / T_earth)

def synodic_period_superior(T_planet, T_earth):
    """
    Calculate synodic period for a superior planet.
    T_planet: orbital period of planet (days)
    T_earth: orbital period of Earth (days)
    Returns: synodic period (days)
    """
    return 1 / (1 / T_earth - 1 / T_planet)

planet_periods = {
    'Mercury': 87.969,
    'Venus': 224.701,
    'Earth': 365.256,
    'Mars': 686.98,
    'Jupiter': 4332.589,
    'Saturn': 10759.22,
    'Uranus': 30685.4,
    'Neptune': 60189,
    'Pluto': 90560
}

T_earth = planet_periods['Earth']

# Calculate synodic periods for all planets except Earth itself
for planet, period in planet_periods.items():
    if planet == 'Earth':
        continue
    if period < T_earth:
        S = synodic_period_inferior(period, T_earth)
        print(f"{planet}: (inferior) synodic period = {S:.2f} days")
    else:
        S = synodic_period_superior(period, T_earth)
        print(f"{planet}: (superior) synodic period = {S:.2f} days")
