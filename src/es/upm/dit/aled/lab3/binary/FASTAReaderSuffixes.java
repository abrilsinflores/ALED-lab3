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
		int hi = this.suffixes.length; // mayor posicc
		boolean found = false;
		int index = 0; // rastrea el caracter actual con el pattern // el FASTA file
		int posSuffix;
		int posInSuffixes = 0;
		int m;

		while ((!found) && ((hi - lo) > 1)) {
			m = (int) Math.floor((lo + hi) / 2); // cálculo punto medio intervalo
			posSuffix = suffixes[m].suffixIndex; // posición del suffix ( q está en la posicc m del array the

			// If the pattern value at 'index' IS NOT THE LAST CHARACTER and matches with
			// the one in the suffix, test the next one
			if (pattern[index] == content[posSuffix + index]) {
				index++;
			}

			// If the pattern value at 'index' IS THE LAST CHARACTER and matches with
			// the one in the suffix, we have found a matching sequence
			if (index == pattern.length && (pattern[index - 1] == content[posSuffix + index - 1])) {
				listOfConcurrences.add(posSuffix);
				found = true;
				posInSuffixes = m;
			}

			// If the pattern value at 'index' comes before the one in the suffix...
			else if (pattern[index] < content[posSuffix + index]) {
				// ...take the left half of the suffix list
				// and restart the index
				hi = m--;
				index = 0;
			}
			// If the pattern value at 'index' comes after the one in the suffix...
			else if (pattern[index] > content[posSuffix + index]) {
				// ...take the right half of the suffix list
				// and restart the index
				lo = m++;
				index = 0;
			}
		}
		if (found) {
			// Now we also check the previous indexes, in case there are more matches
			int indexSubstract = 1;
			while (true) {
				posSuffix = suffixes[posInSuffixes - indexSubstract].suffixIndex;
				boolean isAlsoMatch = true;
				for (int s = 0; s < pattern.length; s++) {
					if (pattern[s] != content[posSuffix + s]) {
						isAlsoMatch = false;
						break;
					}
				}
				if (isAlsoMatch) {
					listOfConcurrences.add(posSuffix);
					indexSubstract++;
				} else
					break;
			}
			// Now we also check the next indexes, in case there are more matches
			int indexAdd = 1;
			while (true) {
				posSuffix = suffixes[posInSuffixes + indexAdd].suffixIndex;
				boolean isAlsoMatch = true;
				for (int s = 0; s < pattern.length; s++) {
					if (pattern[s] != content[posSuffix + s]) {
						isAlsoMatch = false;
						break;
					}
				}
				if (isAlsoMatch) {
					listOfConcurrences.add(posSuffix);
					indexAdd++;
				} else
					break;
			}
		}
		return listOfConcurrences;
	}

	// ACTIVIDAD EXTRA -> implementar searchSNV
	public List<Integer> searchSNV(byte[] pattern) {
		List<Integer> listOfConcurrences = new ArrayList<Integer>();
		int lo = 0;
		int hi = this.suffixes.length - 1;
		boolean found = false;
		int index = 0;

		while ((!found) && !(hi - lo <= 1)) {
			int m = (int) Math.floor(lo + (hi - lo) / 2);
			int posSuffix = suffixes[m].suffixIndex;

			while ((index < pattern.length) && (posSuffix + index < content.length)
					&& (pattern[index] == content[posSuffix + index])) {
				index++;
			}
			// pattern.length-1 porq permitimos un SNV
			if ((index == pattern.length || index == pattern.length - 1) && (posSuffix + index == content.length)) {

				listOfConcurrences.add(posSuffix);

				found = true;

			} else if ((index == pattern.length) || (pattern[index] < content[posSuffix + index])) {
				hi = m--;
				index = 0;
			} else {
				lo = m++;
				index = 0;
			}
		}
		int posSuffix = listOfConcurrences.get(0);
		for (int posNextSuffix = posSuffix + 1; posNextSuffix < this.suffixes.length; posNextSuffix++) {
			index = 0;
			while ((index < pattern.length) && (pattern[index] == this.content[posNextSuffix + index])) {
				index++;
			}
			if (index == pattern.length || index == pattern.length - 1) {
				listOfConcurrences.add(posNextSuffix);
			}
		}
		for (int posBeforeSuffix = posSuffix - 1; posBeforeSuffix > 0; posBeforeSuffix--) {
			index = 0;
			while ((index < pattern.length) && (pattern[index] == this.content[posBeforeSuffix + index])) {
				index++;
			}
			if (index == pattern.length || index == pattern.length - 1) {
				listOfConcurrences.add(posBeforeSuffix);
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
		List<Integer> posiciones = reader.searchSNV(patron);
		System.out.println("Tiempo de búsqueda: " + (System.nanoTime() - t3));
		if (posiciones.size() > 0) {
			for (Integer pos : posiciones)
				System.out.println("Encontrado " + args[1] + " en " + pos);
		} else
			System.out.println("No he encontrado " + args[1] + " en ningún sitio.");
		System.out.println("Tiempo total: " + (System.nanoTime() - t1));
	}
}
