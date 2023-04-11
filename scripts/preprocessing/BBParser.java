import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.io.FileWriter;

/**
 * @author Saurav Kiri
 * @since 2023-04
 */

public class BBParser {
    
    public static void main(String[] args) {
        ArrayList<String> bblines = new ArrayList<>();
        String strLine;

        try (
            FileReader fr = new FileReader(args[0]);
            BufferedReader br = new BufferedReader(fr);
        ) {
            while ((strLine = br.readLine()) != null) {
                bblines.add(strLine);
            }
        } catch (IOException io) {
            System.out.println(io);
        }

        // keep only the lines that start from Input: to Result:, as these are the main BBDuk results
        ArrayList<Integer> inputResIndices = new ArrayList<>();
        for (int i = 0; i < bblines.size(); i++) {
            if (bblines.get(i).startsWith("Input:")) {
                inputResIndices.add(i);
            } else if (bblines.get(i).startsWith("Result:")) {
                inputResIndices.add(i);
            }
        }

        // ensure inputResIndices is sorted in ascending order to obtain the proper subList object
        Collections.sort(inputResIndices);
        List<String> bbsubset;
        try {
            bbsubset = bblines.subList(inputResIndices.get(0), inputResIndices.get(1));
            // write ArrayList to file using FileWriter object
            if (args.length != 2) {
                System.out.println("WARNING: Script requires 2 arguments - BBDuk file and output filename. Subsetted file NOT WRITTEN.");
            } else {
                try (FileWriter fw = new FileWriter(args[1])) {
                    for (String bbstring: bbsubset) {
                        fw.write(bbstring + System.lineSeparator());
                    }
                } catch (IOException io) {
                    System.out.println(io);
                }
            }
        } catch (IndexOutOfBoundsException obe) {
            System.out.println(obe);
            System.exit(1);
        }
    }
}