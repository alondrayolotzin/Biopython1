from Bio import SeqIO

def encontrar_codones(marco, secuencia):
    """
    Encuentra todos los codones en el marco de lectura especificado.
    """
    codones = []
    for i in range(marco, len(secuencia), 3):
        codon = secuencia[i:i+3]
        if len(codon) == 3:
            codones.append(codon)
    return codones

def escribir_codones_a_archivo(marco, secuencia, archivo_salida, hebra):
    """
    Escribe los codones en formato FASTA en el archivo especificado.
    """
    codones = encontrar_codones(marco, secuencia)
    with open(archivo_salida, 'w') as archivo:
        archivo.write(f">Frame_{marco+1}_{hebra}\n")
        archivo.write(" ".join(codones) + "\n")

def procesar_marcos(secuencia, nombre_base_archivo, hebra, secuencia_id):
    """
    Procesa los 6 marcos de lectura y escribe los codones en archivos separados.
    """
    for marco in range(3):
        archivo_salida = f"{nombre_base_archivo}_seq{secuencia_id}_Frame{marco+1}_{hebra}.fa"
        escribir_codones_a_archivo(marco, secuencia, archivo_salida, hebra)

# Cargar las secuencias desde el archivo FASTA
ruta_archivo = '/Users/alondramarquez/Desktop/Biopython1/archivospython/files/seq.nt.fa'  # Cambia esto a la ruta de tu archivo local
secuencias = list(SeqIO.parse(ruta_archivo, "fasta"))

# Procesar cada secuencia en el archivo
for i, registro in enumerate(secuencias):
    secuencia = str(registro.seq)
    secuencia_reversa = secuencia[::-1]

    # Procesar los marcos de lectura para la hebra original
    procesar_marcos(secuencia, "Frame", "directa", i+1)

    # Procesar los marcos de lectura para la hebra inversa
    procesar_marcos(secuencia_reversa, "Frame", "inversa", i+1)