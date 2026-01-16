import pandas as pd
import glob
import os

# Путь к каталогу с результатами Kallisto
base_path = "/home/nicolaedrabcinski/research/metbio/results/viorel_dataset_KALLISTO/quantification"

# Найти все файлы abundance.tsv во всех sample-папках
tsv_files = sorted(glob.glob(os.path.join(base_path, "sample_*", "abundance.tsv")))

# Список для объединённых таблиц
dfs = []

for file_path in tsv_files:
    sample_name = os.path.basename(os.path.dirname(file_path))
    df = pd.read_csv(file_path, sep="\t")
    df["sample"] = sample_name  # добавить имя сэмпла как столбец
    dfs.append(df)

# Объединить все таблицы
combined_df = pd.concat(dfs, ignore_index=True)

# Посмотреть, какие lineages (транскрипты) встречаются
print("\n✅ Найдены колонки:", combined_df.columns.tolist())
print("\n✅ Всего строк:", len(combined_df))
print("\n✅ Примеры данных:")
print(combined_df.head())

# Проверим уникальные lineages
if "target_id" in combined_df.columns:
    print("\n🧬 Уникальные lineages (target_id):")
    print(sorted(combined_df["target_id"].unique()))
else:
    print("\n⚠️ Колонка 'target_id' отсутствует — проверь формат abundance.tsv.")

who_names_mapping = {
    'B.1.1.7': 'Alpha', 'B.1.351': 'Beta', 'P.1': 'Gamma', 'B.1.617.2': 'Delta',
    'B.1.1.529': 'Omicron', 'BA.1': 'Omicron BA.1', 'BA.2': 'Omicron BA.2',
    'BA.3': 'Omicron BA.3', 'BA.4': 'Omicron BA.4', 'BA.5': 'Omicron BA.5',
    'BQ.1': 'Omicron BQ.1', 'XBB': 'Omicron XBB', 'XBB.1.5': 'Omicron XBB.1.5',
    'EG.5': 'Omicron EG.5', 'JN.1': 'Omicron JN.1', 'KP.2': 'Omicron KP.2', 'KP.3': 'Omicron KP.3',
    'B.1.429': 'Epsilon', 'B.1.427': 'Epsilon', 'P.2': 'Zeta',
    'B.1.525': 'Eta', 'B.1.526': 'Iota', 'B.1.617.1': 'Kappa',
    'C.37': 'Lambda', 'B.1.621': 'Mu', 'P.3': 'Theta'
}

present = {'B.1.1.7','B.1.351','B.1.429','B.1.525','B.1.526','B.1.617.1','B.1.617.2','C.37','P.1','P.2','P.3'}

missing = set(who_names_mapping.keys()) - present
print("Отсутствующие lineages:", missing)
