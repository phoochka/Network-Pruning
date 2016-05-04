package pruning;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public class Edge extends Pair {


    private int node1 = (int) getA();
    private int node2 = (int) getB();

    public Edge(int node1, int node2){
        super(node1, node2);
    }

    public int getNode1() {
        return node1;
    }

    public int getNode2() {
        return node2;
    }
    public boolean containsNode(int node) {
        return node == node1 || node == node2;
    }

    public boolean isConnected(Set<Edge> searchSet) {
        for (Edge e: searchSet) {
            if (e.containsNode(this.node1) || e.containsNode(this.node2)) return true;
        }
        return false;
    }

    public Set<Edge> getAnyConnected(Set<Edge> searchSet) {
        return searchSet.parallelStream()
                .filter(edge -> edge.containsAny(this))
                .collect(Collectors.toSet());
    }

    public boolean listHasEdge(Set<Integer> searchSet) {
        return searchSet.contains(node1) && searchSet.contains(node2);
    }

    public boolean listHasPartialEdge(Set<Integer> searchSet) {
        if (listHasEdge(searchSet)) return false;
        else return searchSet.contains(node1) || searchSet.contains(node2);
    }

    public boolean containsAny(Edge searchEdge) {
        return searchEdge.containsNode(this.node1) || searchEdge.containsNode(this.node2);
    }

    public boolean containsAny(Set<Edge> searchSet) {
        return searchSet.parallelStream().anyMatch(edge -> edge.equals(this));
    }

    public int getOtherNode(int node) {
        if (node == node1) return node2;
        else if (node == node2) return node1;
        else return 0;
    }

    public Set<Integer> toSet() {
        Set<Integer> result = new HashSet<Integer>(2);
        result.add(node1);
        result.add(node2);
        return result;
    }

    @Override
    public String toString() {
        return (node1+","+node2);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Edge edgePair = (Edge) o;

        if (node1 == edgePair.node1 && node2 == edgePair.node2) return true;
        else return (node2 == edgePair.node1 && node1 == edgePair.node2);
    }

    @Override
    public int hashCode() {
        int result = node1 * node2;
        result = result + (31 * (node1 + node2));
        return result;
    }
}
