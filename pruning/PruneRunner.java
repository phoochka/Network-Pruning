package pruning;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;

public class PruneRunner {

    public PruneRunner() {
        String filename = "UBD080UBC040Gv001000Ge20Gg0100Bv0010Bd0005Bg01000400_10.txt";
        Prune prune = new Prune(filename);

        double[] kList = {1.0,1.2,1.4,1.5,1.6,1.8,2};

        //outputKFactor(prune, kList);

        for (double k :kList) {
            prune.traverseGraph(prune.getConnectedComponents(k));
        }

    }

    public void outputKFactor(Prune prune, double[] kList) {
        try {
            PrintWriter writer = new PrintWriter("thresholding_result.txt");
            for (double k : kList) {
                writer.println("Size of list, Average size of connected nodes, for k = " + k);
                for (List<Set<Edge>> l : prune.getConnectedComponents(k).values()) {
                    writer.println(l.size()+","+l.stream().mapToInt(set -> set.size()).average().getAsDouble());
                }
            }
            writer.close();
        } catch (FileNotFoundException f) {
            System.out.println("Could not write");
        }
    }

    public static void main(String[] qwerty) {



        new PruneRunner();
    }
}
