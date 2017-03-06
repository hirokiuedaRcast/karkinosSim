package jp.ac.utokyo.rcast.karkinosSim.util;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Linked hash map that limits number of added entries. If exceeded the eldest
 * element is removed.
 * 
 * @param <K>
 *            key type
 * @param <V>
 *            value type
 */
public class LimitedLinkedHashMap<K, V> extends LinkedHashMap<K, V> {

	private final int maxSize;

	/**
	 * Constructor.
	 * 
	 * @param maxSize
	 *            maximum number of elements given hash map can hold.
	 */
	public LimitedLinkedHashMap(int maxSize) {
		this.maxSize = maxSize;
	}

	@Override
	protected boolean removeEldestEntry(java.util.Map.Entry<K, V> eldest) {
		return maxSize != 0 && size() > maxSize;
	}


}
