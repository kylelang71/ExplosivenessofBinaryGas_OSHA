import pandas as pd
import PR_EOS
import os

# LOAD DATA
path = os.getcwd()
raw_properties = pd.read_csv(f'{path}/Physical_Properties.csv')
explosive_limits = pd.read_csv(f'{path}/Explosive_Limits.csv')
product_recipe = pd.read_csv(f'{path}/Glass_Cleaner.csv')

# PULL / STRUCTURE PROPERTIES
recipe = product_recipe[['NAME', 'X']]
df = pd.merge(recipe, raw_properties, how='left')
df = pd.merge(df, explosive_limits, how='left')

print(df)
print("")

# CALCULATE MOLAR VAPOR COMPOSITION W/ PENG-ROBINSON EOS
y = PR_EOS.molar_vapor_composition_STP(df.at[0, 'omega'], df.at[0, 'Tc(K)'], df.at[0, 'Pc(MPa)'], df.at[0, 'X'],
                                       df.at[1, 'omega'], df.at[1, 'Tc(K)'], df.at[1, 'Pc(MPa)'], df.at[1, 'X'])
y = {key: round(value, 5) for key, value in y.items()}
y1 = y.get('Comp1_MoleFraction')
y2 = y.get('Comp2_MoleFraction')

print('Molar Vapor Composition of Ingredient 1:')
print(y1)
print('Molar Vapor Composition of Ingredient 2:')
print(y2)
print()

for row in df:
    if df.at[0, 'LEL'] <= (100 * y1) and df.at[0, 'UEL'] >= (100 * y1):
        print("Component 1's concentration in the product's vapor is within it's flammability limits "
              "under standard conditions.  It may be flammable without proper ventilation.")
        break
    elif df.at[1, 'LEL'] <= (100 * y2) and df.at[1, 'UEL'] >= (100 * y2):
        print("Component 2's concentration in the product's vapor is within it's flammability limit "
              "under standard conditions. It may be flammable without proper ventilation.")
        break
    else:
        print("The product's vapor likely is not flammable.")
        break
