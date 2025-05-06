import os
import subprocess
import argparse
import pandas as pd

def check_fasta_not_empty(fasta_file):
    """
    Comprueba si el archivo FASTA no está vacío y contiene secuencias.
    """
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                return True
    return False

def run_hmmscan(hmm_file, fasta_file, output_file,log_file):
    """
    Ejecuta hmmscan con los argumentos especificados.
    """
    subprocess.run([
        #"/beegfs/data/soft_legacy/hmmer-3.1b2/bin/hmmscan",
        "hmmscan",
        "-o", log_file,        
        "--tblout", output_file,
        hmm_file,
        fasta_file
    ], check=True)

def extract_second_word(target_name):
    """
    Extrae la segunda palabra del nombre del target para simplificar.
    """
    parts = target_name.split("_")
    second_word = parts[1] if len(parts) > 1 else parts[0]
    return second_word

def parse_hmmscan_output(output_file):
    """
    Analiza la salida de hmmscan y genera una tabla con cuatro columnas:
    - Target name
    - Query name
    - Best score
    - Score ratio
    """
    results = {}
    with open(output_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            try:
                query_name = parts[2]
                #target_name = extract_second_word(parts[0])
                target_name = parts[0]
                full_sequence_score = float(parts[5])

                if query_name not in results:
                    results[query_name] = []
                results[query_name].append((target_name, full_sequence_score))
            except (IndexError, ValueError):
                continue

    formatted_results = []
    for query, hits in results.items():
        hits.sort(key=lambda x: x[1], reverse=True)
        best_hit = hits[0] if hits else (None, 0)
        second_hit = hits[1] if len(hits) > 1 else (None, 0)

        target_name = best_hit[0] if best_hit[0] else "No hit"
        best_score = best_hit[1] if best_hit[1] else 0
        score_ratio = (
            best_hit[1] / second_hit[1] if second_hit[1] > 0 else "no second hit"
        )
        formatted_results.append((target_name, query, best_score, score_ratio))

    return formatted_results

def write_output_table(results, output_table_file):
    """
    Escribe la tabla formateada en un archivo de salida.
    """
    with open(output_table_file, "w") as f:
        f.write("Target Name\tQuery Name\tBest Score\tScore Ratio\n")
        for row in results:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")

def update_csv_with_results(csv_file, formatted_results):
    """
    Crea un nuevo archivo CSV actualizado con los resultados obtenidos de hmmscan.
    """
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file, sep=";", index_col=0)
        # Traitement specifique genewise (crade)
        if "Stop/Shift Positions" in df:
            df = pd.read_csv(csv_file, sep=";", index_col=0,dtype={"Stop/Shift Positions": str})
        if df.empty:
            return

        # Eliminar filas duplicadas basadas en la columna SeqID, conservando solo la primera
        df = df[~df.duplicated(subset='SeqID', keep='first')]
        for target_name, query, best_score, score_ratio in formatted_results:
            #match_row = df[df['SeqID'].str.contains(query, na=False)]
            match_row = df[df['SeqID'].str.fullmatch(query, na=False)]
            if not match_row.empty:
                index = match_row.index[0]
                df.at[index, 'Best Match'] = target_name
                df.at[index, 'Bit Score'] = best_score
                df.at[index, 'Score ratio'] = score_ratio
            else:
                print("DEBUG SIMON " + target_name +" "+ query)
        # Crear el nuevo archivo con un sufijo "_curated"
        curated_csv_file = os.path.splitext(csv_file)[0] + "_curated.csv"
        df.to_csv(curated_csv_file, sep=";", encoding="utf-8")

def main():
    parser = argparse.ArgumentParser(description="Ejecuta hmmscan y genera una tabla con score ratios y best scores.")
    parser.add_argument("--hmm_db", required=True, help="Archivo de base de datos HMM")
    parser.add_argument("--input_fasta", required=True, help="Archivo de entrada en formato FASTA")
    parser.add_argument("--output_name", required=True, help="Nombre del archivo de salida para hmmscan")
    parser.add_argument("--output_table", required=True, help="Nombre del archivo de tabla formateada")
    parser.add_argument("--csv_file", help="Archivo CSV con los datos a actualizar")

    args = parser.parse_args()

    if not check_fasta_not_empty(args.input_fasta):
        # Si el archivo FASTA está vacío, generar los archivos vacíos
        open(args.output_name, 'w').close()  # Crear archivo de salida vacío para hmmscan
        write_output_table([], args.output_table)  # Crear tabla vacía

        if args.csv_file:
            # Si se pasa un archivo CSV, crear un archivo CSV vacío
            df = pd.read_csv(args.csv_file, sep=";", index_col=0)
            df.to_csv(os.path.splitext(args.csv_file)[0] + "_curated.csv", sep=";", encoding="utf-8")
        return

    # Si el FASTA no está vacío, proceder con el análisis normal
    logfile = os.path.splitext(args.output_name)[0] + ".log"
    print("HMMSCAN HMMFILE = " + args.hmm_db + " FASTA = " + args.input_fasta)
    run_hmmscan(args.hmm_db, args.input_fasta, args.output_name,logfile)
    results = parse_hmmscan_output(args.output_name)
    write_output_table(results, args.output_table)

    if args.csv_file:
        update_csv_with_results(args.csv_file, results)

if __name__ == "__main__":
    main()

