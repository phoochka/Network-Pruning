package pruning.util;

import pruning.Edge;

import java.util.HashMap;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by gaurav on 4/24/16.
 */
public class GraphStats {

    private HashMap<Edge, Integer> squashedGraph;
    private int size;
    private double mean;
    private double stdDev;

    private void getStats() {
        size = squashedGraph.values().size();

        // Finding mean
        // int sum = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).sum();
        mean = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).average().getAsDouble();

        // System.out.println("Mean is: "+mean);

        // Finding standard deviation
        double temp = 0;
        for (int i : squashedGraph.values()) {
            temp += (mean - i) * (mean - i);
        }
        stdDev = Math.sqrt(temp / size);
        // System.out.println("stddev is: "+stdDev);
    }

    public Set<Edge> thresholdEdgesByStd(double k) {
        double thresh = mean + (k * stdDev);
        return squashedGraph.keySet().parallelStream()
            .filter(e -> squashedGraph.get(e) > thresh).collect(Collectors.toSet());
    }


}
