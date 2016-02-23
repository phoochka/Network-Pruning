package pruning;

import java.util.Set;

public class Result {

    double normalizedConductance;
    double rawConductance;
    Set<Integer> nodes;
    int startTime;
    int endTime;

    public Result(double normalizedConductance, double rawConductance, Set<Integer> nodes, int startTime, int endTime) {
        this.normalizedConductance = normalizedConductance;
        this.rawConductance = rawConductance;
        this.nodes = nodes;
        this.startTime = startTime;
        this.endTime = endTime;
    }

    @Override
    public String toString() {
        return "\nNormalized Conductance:"+normalizedConductance
                +"Raw Conductance: "+rawConductance
                +"\nStart Time: "+startTime+" End Time: "+endTime
                +"\nNodes: "+nodes;
    }
}
