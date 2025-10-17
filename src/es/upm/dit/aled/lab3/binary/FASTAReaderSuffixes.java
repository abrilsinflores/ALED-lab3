package es.upm.dit.aled.lab3.binary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import es.upm.dit.aled.lab3.FASTAReader;

/**
 * Reads a FASTA file containing genetic information and allows for the search
 * of specific patterns within these data. The information is stored as an array
 * of bytes that contain nucleotides in the FASTA format. Since this array is
 * usually created before knowing how many characters in the origin FASTA file
 * are valid, an int indicating how many bytes of the array are valid is also
 * stored. All valid characters will be at the beginning of the array.
 * 
 * This extension of the FASTAReader uses a sorted dictionary of suffixes to
 * allow for the implementation of binary search.
 * 
 * @author mmiguel, rgarciacarmona
 *
 */
public class FASTAReaderSuffixes extends FASTAReader { // tiene tous argumentos FASTAReader + array sufijos
	protected Suffix[] suffixes;

	/**
	 * Creates a new FASTAReader from a FASTA file.
	 * 
	 * At the end of the constructor, the data is sorted through an array of
	 * suffixes.
	 * 
	 * @param fileName The name of the FASTA file.
	 */
	public FASTAReaderSuffixes(String fileName) {
		// Calls the parent constructor
		super(fileName);
		this.suffixes = new Suffix[validBytes];
		for (int i = 0; i < validBytes; i++)
			suffixes[i] = new Suffix(i); // crea tantos sufijos como Bytes tenga el archivo
		// Sorts the data
		sort();
	}

	/*
	 * Helper method that creates an array of integers that contains the positions
	 * of all suffixes, sorted alphabetically by the suffix.
	 */
	private void sort() {
		// Instantiate the external SuffixComparator, passing 'this' (the reader)
		// so it can access the content and validBytes fields.
		SuffixComparator suffixComparator = new SuffixComparator(this);
		// Use the external Comparator for sorting.
		Arrays.sort(this.suffixes, suffixComparator);
	}

	/**
	 * Prints a list of all the suffixes and their position in the data array.
	 */
	public void printSuffixes() {
		System.out.println("-------------------------------------------------------------------------");
		System.out.println("Index | Sequence");
		System.out.println("-------------------------------------------------------------------------");
		for (int i = 0; i < suffixes.length; i++) {
			int index = suffixes[i].suffixIndex;
			String ith = "\"" + new String(content, index, Math.min(50, validBytes - index)) + "\"";
			System.out.printf("  %3d | %s\n", index, ith);
		}
		System.out.println("-------------------------------------------------------------------------");
	}

	/**
	 * Implements a binary search to look for the provided pattern in the data
	 * array. Returns a List of Integers that point to the initial positions of all
	 * the occurrences of the pattern in the data.
	 * 
	 * @param pattern The pattern to be found.
	 * @return All the positions of the first character of every occurrence of the
	 *         pattern in the data.
	 */
	@Override
	public List<Integer> search(byte[] pattern) { // TODO

		List<Integer> listOfConcurrences = new ArrayList<Integer>();
		// límites de la búsqueda binaria
		int lo = 0;
		int hi = this.suffixes.length - 1; // mayor posicc (-1 xq si mide 8 tengo hasta la posicc 7!!!)
		boolean found = false;
		int index = 0; // rastrea el caracter actual con el pattern

		while ((!found) && !(hi - lo <= 1)) {
			int m = (int) Math.floor(lo + (hi - lo) / 2); // cálculo punto medio intervalo
			int posSuffix = suffixes[m].suffixIndex; // posición del suffix ( q está en la posicc m del array the
														// suffixs)en el FASTA file
			// comparacc de cada carácter con el pattern mientras que:
			// 1) el índice no supere el tamaño del patrón
			// 2) los bytes del suffix q miremos no se salgan del archivo
			// 3) los bytes del suffix y el archivo coincidan
			while ((index < pattern.length) && (posSuffix + index < content.length)
					&& (pattern[index] == content[posSuffix + index])) {
				index++; // aumentamos el índice
			}

			// añadimos la posición del sufijo a la lista SOLAMENTE si
			// coinciden el nm d bytes del suffix iguales al del pattern
			// y si NO HAY MÁS LETRAS DETRÁS
			if ((index == pattern.length) && (posSuffix + index == content.length)) {
				listOfConcurrences.add(posSuffix);
				found = true;
			}
			// si el suffix no coincide, MODIFICAMOS LÍMITES INTERVALOS
			// si la letra del sufijo es MAYOR q la del pattern, CONCIDICIONES:
			// 1) q el índice sea igual d largo (xq entonces se ha encontrado un patrón pero
			// DESPUÉS HAY MÁS LETRAS)
			// OOOOO 2) q la sig letra del sufijo sea MAYOR q la del patrón
			// OJO:can comparar letras porq es código ASCII, está ya ordenado alfabéticamm
			// según Byte que es
			else if ((index == pattern.length) || (pattern[index] < content[posSuffix + index])) {
				hi = m--;
				index = 0;
			} else {
				lo = m++;
				index = 0;
			}
		}

		return listOfConcurrences;
	}

	public static void main(String[] args) {
		long t1 = System.nanoTime();
		FASTAReaderSuffixes reader = new FASTAReaderSuffixes(args[0]);
		if (args.length == 1)
			return;
		byte[] patron = args[1].getBytes();
		System.out.println("Tiempo de apertura de fichero: " + (System.nanoTime() - t1));
		long t2 = System.nanoTime();
		System.out.println("Tiempo de ordenación: " + (System.nanoTime() - t2));
		reader.printSuffixes();
		long t3 = System.nanoTime();
		List<Integer> posiciones = reader.search(patron);
		System.out.println("Tiempo de búsqueda: " + (System.nanoTime() - t3));
		if (posiciones.size() > 0) {
			for (Integer pos : posiciones)
				System.out.println("Encontrado " + args[1] + " en " + pos);
		} else
			System.out.println("No he encontrado " + args[1] + " en ningún sitio.");
		System.out.println("Tiempo total: " + (System.nanoTime() - t1));
	}
}
