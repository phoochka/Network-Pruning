package lkml_reply;


import pruning.Edge;
import pruning.Prune;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by gaurav on 5/30/16.
 */
public class ProcessMailingListData {

    final static long EPOCH_MONTH = 2629743;

    private HashMap<Edge, Integer> squashedGraph;

    private int size;
    private double mean;
    private double stdDev;

    public ProcessMailingListData() {
        TreeMap<Integer, HashMap<Edge, Integer>> timeGraph = processFile(new File("lkml_reply/out.lkml-reply"));
//        connectingGraph(timeGraph);
//        checkConnectivity(timeGraph);

        System.out.println("Average no. of edges per time period: "+timeGraph.values().stream().mapToInt(map-> map.size()).average());

        getStats();
        TreeMap<Integer, List<Set<Edge>>> prunedMap = getConnectedComponents(100, timeGraph);

        findIntersection(prunedMap);


//        writeToFile(timeGraph);
    }

    private void findIntersection(TreeMap<Integer, List<Set<Edge>>> prunedMap) {

        Set<Edge> common;

        for (int range = 1; range < prunedMap.size() - 1; range++) {
            for(int n = prunedMap.firstKey(); n < (prunedMap.lastKey() - range); n++){
                int time = n;
                int end = n + range;

                common = new HashSet<>();
                prunedMap.get(time).forEach(common::addAll);

                for (;time<end;time++) {
                    Set<Edge> edgesInThisTime = new HashSet<>();
                    prunedMap.get(time).forEach(set -> edgesInThisTime.addAll(set));

//                if (time == prunedMap.lastKey()) System.out.println("Size at time " + time + ": " + common.size());
                        common.retainAll(edgesInThisTime);
//                if (time == prunedMap.lastKey()) System.out.println("After prune  " + time + ": " + common.size());
                }

                System.out.println("From "+(n-prunedMap.firstKey())+" to "+(end-prunedMap.firstKey())+" intersections: "+common.size());
            }
        }
    }

    private TreeMap<Integer, List<Set<Edge>>> getConnectedComponents(double k, TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {

        TreeMap<Integer, List<Set<Edge>>> result = new TreeMap<>();
        for (int time : timeGraph.keySet()) {
//            System.out.println("processing time period "+time);

            HashMap<Edge, Integer> edgeMap = timeGraph.get(time);

            Set<Edge> thresholdedEdges = thresholdEdgesByPer(k, edgeMap);

            List<Set<Edge>> timeList = new ArrayList<Set<Edge>>();
            Set<Edge> connectedSet;

            Set<Integer> nodesToSearch = new HashSet<>();
            edgeMap.keySet().stream()
//                    .filter(edge -> thresholdedEdges.contains(edge))
                    .forEach(edge -> nodesToSearch.addAll(edge.toSet()));

            while (!nodesToSearch.isEmpty()) {
                Stack<Integer> stack = new Stack<>();
//                System.out.println("nodes to search: "+nodesToSearch+" "+nodesToSearch.isEmpty()+" "+nodesToSearch.size());
//                if(nodesToSearch.size() == 3) nodesToSearch.forEach(node -> System.out.println("n "+node));
                stack.push(nodesToSearch.stream().findFirst().get());
                connectedSet = new HashSet<Edge>();
                while (!stack.isEmpty()) {
                    Integer n1 = stack.pop();
                    nodesToSearch.remove(n1);
                    Set<Edge> edgesToSearch = edgeMap.keySet().stream()
                            .filter(thresholdedEdges::contains).collect(Collectors.toSet());

                    Edge edge = new Edge(n1, n1);
                    for (Edge connectedEdge : edge.getAnyConnected(edgesToSearch)) {
                        if (!connectedSet.contains(connectedEdge)) {
                            connectedSet.add(connectedEdge);
                            stack.push(n1);
                            nodesToSearch.remove(n1);
                        }
                    }
                }
                timeList.add(connectedSet);
            }
            result.put(time, timeList);
        }
        return result;
    }

    public Set<Edge> thresholdEdgesByPer(double k, HashMap<Edge, Integer> edgeWeights) {
        int edgeListSize = (int) (edgeWeights.size() * (k / 100));
        Comparator<Map.Entry<Edge, Integer>> sortEdgesByWeight = (e1, e2) -> e1.getValue().compareTo(e2.getValue());
        return squashedGraph.entrySet().stream()
                .sorted(sortEdgesByWeight.reversed())
                .limit(edgeListSize)
                .map(e -> e.getKey())
                .collect(Collectors.toSet());
    }

    private void checkConnectivity(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {
        for (int time : timeGraph.keySet()) {
            HashMap<Edge, Integer> edgeMap = timeGraph.get(time);
            Set<Edge> connectedSet;

            List<Set<Edge>> timeList = new ArrayList<Set<Edge>>();
            Set<Edge> nodesToSearch = new HashSet<Edge>();
            nodesToSearch.addAll(edgeMap.keySet().stream()
                    .collect(Collectors.toSet()));

            while (!nodesToSearch.isEmpty()) {
                Stack<Edge> stack = new Stack<Edge>();
                stack.push(nodesToSearch.stream().findFirst().get());
                connectedSet = new HashSet<Edge>();
                while (!stack.isEmpty()) {
                    Edge edge = stack.pop();
                    nodesToSearch.remove(edge);
                    Set<Edge> edgesToSearch = edgeMap.keySet().stream()
                            .collect(Collectors.toSet());

                    for (Edge connectedEdge : edge.getAnyConnected(edgesToSearch)) {
                        if (!connectedSet.contains(connectedEdge)) {
                            connectedSet.add(connectedEdge);
                            stack.push(connectedEdge);
                            nodesToSearch.remove(connectedEdge);
                        }
                    }
                }
                timeList.add(connectedSet);
            }

            System.out.println("Time list size: " + timeList.size());
        }
    }

    private void writeToFile(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {
        try {
            PrintWriter writer = new PrintWriter("mailing_list_data.txt");
            writer.println("Node1,Node2,Time,Weight");
            for (int time : timeGraph.keySet()) {
                for (Edge edge : timeGraph.get(time).keySet()) {
                    writer.println(edge.getNode1() + "," + edge.getNode2() + "," + time + "," + timeGraph.get(time).get(edge));
                }
            }
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void connectingGraph(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {

        Set<Edge> allEdges = new HashSet<>();
        timeGraph.values().parallelStream().map(map -> map.keySet()).forEach(e -> allEdges.addAll(e));

        for (int time : timeGraph.keySet()) {
            for (Edge edge : allEdges) {
                if (!timeGraph.get(time).keySet().contains(edge)) timeGraph.get(time).put(edge, 1);
            }
        }
        System.out.println("Finished connecting graph");
    }

    private void getStats() {
        size = squashedGraph.values().size();

        // Finding mean
        // int sum = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).sum();
        mean = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).average().getAsDouble();

         System.out.println("Mean is: "+mean);

        // Finding standard deviation
        double temp = 0;
        for (int i : squashedGraph.values()) {
            temp += (mean - i) * (mean - i);
        }
        stdDev = Math.sqrt(temp / size);
         System.out.println("stddev is: "+stdDev);
    }

    public Set<Edge> thresholdEdgesByStd(double k) {
        double thresh = mean + (k * stdDev);
        return squashedGraph.keySet().parallelStream()
                .filter(e -> squashedGraph.get(e) > thresh).collect(Collectors.toSet());
    }

    private TreeMap<Integer, HashMap<Edge, Integer>> processFile(File f) {

        TreeMap<Integer, HashMap<Edge, Integer>> timeGraph = new TreeMap<>();

        squashedGraph = new HashMap<>();
        String[] line;
        try {
            Scanner rawScanner = new Scanner(f);
            rawScanner.nextLine(); // Skipping headers
            while (rawScanner.hasNextLine()) {
                line = rawScanner.nextLine().split("\t");

                int node1 = Integer.parseInt(line[0]) - 1;
                int node2 = Integer.parseInt(line[1]) - 1;

                if (node1 != node2) {

                    Edge edge = new Edge(node1, node2);

                    Integer weight = Integer.parseInt(line[2]);

                    squashedGraph.compute(edge, (k, v) -> v == null ? 1 : v + 1);

                    long epochTime = Long.parseLong(line[3]);
                    int time = Math.toIntExact(epochTime / (EPOCH_MONTH * 2));


//                    System.out.println("Processing edge " + edge.getNode1()
//                        + " " + edge.getNode2() + " at " + time);

                    if (timeGraph.containsKey(time)) {
                        if (timeGraph.get(time).containsKey(edge)) {
                            Integer oldWeight = timeGraph.get(time).get(edge);
                            timeGraph.get(time).replace(edge, oldWeight, weight + oldWeight);
                        } else timeGraph.get(time).put(new Edge(node1, node2), weight);
                    } else {
                        HashMap<Edge, Integer> edgeWeights = new HashMap<>();
                        edgeWeights.put(edge, weight);
                        timeGraph.put(time, edgeWeights);
                    }
                }
            }
            rawScanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        System.out.println("Finished creating data structure total time periods: " + timeGraph.keySet().size());

        return timeGraph;
    }

    public static void main(String[] sdff) {
        new ProcessMailingListData();
    }
}

