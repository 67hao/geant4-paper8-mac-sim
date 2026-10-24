"""
Script don gian: loc cac cot N (so su kien theo tung loai tuong tac)
tu file CSV ket qua mo phong, roi xuat ra file Excel theo bo cuc:

            S0                              S1                        ...
Energy_keV  N_compton N_photo N_pair ...    N_compton N_photo N_pair ...
15
20
...

- Moi vat lieu (S0, S1, S2, ...) la MOT khoi cot rieng, nam canh nhau.
- Hang la nang luong (Energy_keV), lay theo thu tu tang dan.
- Cac cot N duoc lay tu file dau vao: N_compton, N_photo, N_pair,
  N_rayleigh, N_other, N_total (co the sua trong bien N_COLS ben duoi).

Cach dung:
    python loc_N_tuong_tac.py input.csv output.xlsx
"""

import sys
import pandas as pd

# Cac cot "N" (so su kien) muon lay ra - sua lai neu can them/bot
N_COLS = ["N_compton", "N_photo", "N_pair", "N_rayleigh", "N_other", "N_total"]


def build_table(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    # Sap xep de nang luong tang dan, vat lieu theo thu tu xuat hien
    materials = sorted(df["Material"].unique(), key=lambda x: (len(x), x))

    # Bang nang luong duy nhat lam truc hang
    energies = sorted(df["Energy_keV"].unique())
    result = pd.DataFrame({"Energy_keV": energies})

    for mat in materials:
        sub = df[df["Material"] == mat].set_index("Energy_keV")
        for col in N_COLS:
            # Ten cot ket qua: <Material>_<ten cot N>, vi du S1_N_compton
            result[f"{mat}_{col}"] = result["Energy_keV"].map(sub[col])

    return result, materials


def save_excel(result: pd.DataFrame, materials, out_path: str):
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        result.to_excel(writer, sheet_name="N_theo_tuong_tac", index=False, startrow=1)

        ws = writer.sheets["N_theo_tuong_tac"]
        n_types = len(N_COLS)

        # Dong 1: ten vat lieu (S0, S1, ...) gop o tren tung khoi cot
        col = 2  # cot B (cot A la Energy_keV)
        for mat in materials:
            ws.cell(row=1, column=col, value=mat)
            ws.merge_cells(start_row=1, start_column=col, end_row=1, end_column=col + n_types - 1)
            col += n_types

        ws.cell(row=1, column=1, value="")  # o goc tren ben trai de trong
        ws.cell(row=2, column=1, value="Energy_keV")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Cach dung: python loc_N_tuong_tac.py input.csv output.xlsx")
        sys.exit(1)

    csv_path, out_path = sys.argv[1], sys.argv[2]
    result, materials = build_table(csv_path)
    save_excel(result, materials, out_path)
    print(f"Da xuat file: {out_path}")
