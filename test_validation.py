import sys

def validate(Temp, nCo, nNi, V, dt):
    if Temp < 0:
        return False, "Invalid Temp"
    if nCo < 0 or nNi < 0 or (nCo + nNi) <= 0:
        return False, "Invalid Atoms"
    if V <= 0:
        return False, "Invalid V"
    if dt <= 0:
        return False, "Invalid dt"
    return True, "Valid"

tests = [
    (273.15, 5, 5, 109.5, 1e-15, True), # Normal
    (-10, 5, 5, 109.5, 1e-15, False), # Temp negative
    (273.15, -1, 5, 109.5, 1e-15, False), # nCo negative
    (273.15, 0, 0, 109.5, 1e-15, False), # sum = 0
    (273.15, 5, 5, 0, 1e-15, False), # V = 0
    (273.15, 5, 5, -10, 1e-15, False), # V negative
    (273.15, 5, 5, 109.5, 0, False), # dt = 0
    (273.15, 5, 5, 109.5, -1e-15, False), # dt negative
]

all_passed = True
for Temp, nCo, nNi, V, dt, expected in tests:
    result, msg = validate(Temp, nCo, nNi, V, dt)
    if result != expected:
        print(f"Test failed! {Temp, nCo, nNi, V, dt} Expected {expected}, got {result}")
        all_passed = False

if all_passed:
    print("All validation logic tests passed.")
else:
    sys.exit(1)
