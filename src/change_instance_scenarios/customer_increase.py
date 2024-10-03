import pandas as pd

# Lire le fichier txt dans un DataFrame
file_path = 'ton_fichier.txt'
df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt', delim_whitespace=True)

# Demander la valeur de k
# k = float(input("Entrez le pourcentage d'augmentation (par exemple, 10 pour 10%): "))
k = 4

# Calculer la somme des poids avant l'augmentation
sum_before = df['weight'].sum()
print(f"Somme avant augmentation: {sum_before}")

# Augmenter chaque poids de k%
df['weight'] = df['weight'] * (1 + k / 100)

# Calculer la somme des poids après l'augmentation
sum_after = df['weight'].sum()
print(f"Somme après augmentation: {sum_after}")

# Sauvegarder la nouvelle table dans un fichier txt
output_file = 'nouveau_fichier.txt'
df.to_csv(output_file, sep=' ', index=False)

print(f"Table modifiée sauvegardée dans {output_file}")