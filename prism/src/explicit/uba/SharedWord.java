package explicit.uba;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class SharedWord<Letter>
{
	private SharedWord<Letter> prev = null;
	private Letter letter = null;
	private int length = 0;

	/** empty word constructor */
	public SharedWord()
	{
	}

	public SharedWord<Letter> append(Letter letter) {
		SharedWord<Letter> result = new SharedWord<Letter>();
		result.prev = this;
		result.letter = letter;
		result.length = length+1;
		return result;
	}
	
	public int size()
	{
		return length;
	}

	public List<Letter> getWord()
	{
		ArrayList<Letter> list = new ArrayList<Letter>();

		SharedWord<Letter> current = this;
		while (current != null && current.letter != null) {
			list.add(current.letter);
			current = current.prev;
		}

		Collections.reverse(list);
		return list;
	}
	
	public String toString()
	{
		if(length==0) {
			return "ε";
		}
		List<Letter> word = getWord();
		String result = "";
		Iterator<Letter> it = word.iterator();
		result += it.next().toString();
		while(it.hasNext()) {
			result += ";" + it.next().toString();
		}
		return result;
	}
}
