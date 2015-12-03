package pruning;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

class Prune {

    private TreeMap<Integer, HashMap<Edge, Integer>> timeGraph;
    private ConcurrentHashMap<Edge, Integer> squashedGraph;
    private static final int MAX_SEARCH_TIME = 10;

    private int size;
    private double mean;
    private double stdDev;

    private double bestConductance;
    private Set<Integer> bestNodes;
    private int bestStartTime;
    private int bestEndTime;

    private int zeroCount; // TODO: find out why we are getting zero internal weight

    public void resetBestConductance() {
        bestConductance = 0;
    }

    public Prune(String filename) {
        readGraph(filename);
        //printGraph();
        //System.out.println("All nodes at time 82: " + timeGraph.get(82).keySet());
        //System.out.println("All neighbours of 21 at time 82: " + getNeighbours(82, 21));
    }

    private void readGraph(String filename) {
        timeGraph = new TreeMap<Integer, HashMap<Edge, Integer>>();
        squashedGraph = new ConcurrentHashMap<Edge, Integer>();
        String basepath = "output";
        String line;
        try {
            Scanner timeGraphScanner = new Scanner(new File(basepath + "/graph/" + filename));
            timeGraphScanner.nextLine(); // Skip heading line

            while (timeGraphScanner.hasNextLine()) {
                line = timeGraphScanner.nextLine();
                String s[] = line.split(",");
                if (s.length != 3) {
                    System.out.println("Input file in not in expected format");
                    timeGraphScanner.close();
                    throw new NoSuchElementException();
                }
                Integer node1 = Integer.parseInt(s[0]);
                Integer node2 = Integer.parseInt(s[1]);
                Double t = Double.parseDouble(s[2]);
                int floorInt = new Double(Math.floor(t)).intValue();
                Edge nodePair = new Edge(node1, node2);

                squashedGraph.compute(nodePair, (k,v) -> v ==null ? 1 : v + 1);

                HashMap<Edge, Integer> edgeMap;
                if (timeGraph.containsKey(floorInt)) {
                    edgeMap = timeGraph.get(floorInt);
                    edgeMap.compute(nodePair, (k, v) -> v == null ? 1 : v + 1);
                    /**
                     if (edgeMap.containsKey(nodePair)) {
                     edgeMap.put(nodePair, edgeMap.get(nodePair) + 1);
                     } else  edgeMap.put(nodePair, 1);
                     timeGraph.put(floorInt, edgeMap);
                     **/

                } else {
                    edgeMap = new HashMap<Edge, Integer>();
                    edgeMap.put(nodePair, new Integer(1));
                    timeGraph.put(floorInt, edgeMap);
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        getStats();
    }

    private void getStats() {

        size = squashedGraph.values().size();

        System.out.println("Key size: "+squashedGraph.keySet().size()+" value size: "+squashedGraph.values().size());

        // Finding mean
        int sum = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).sum();
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


    public Set<Edge> getThresholdedEdges(double k) {
        double thresh = mean + (k * stdDev);
        return squashedGraph.keySet().parallelStream()
                .filter(e -> squashedGraph.get(e) > thresh).collect(Collectors.toSet());
    }

    public void traverseGraph(TreeMap<Integer, List<Set<Edge>>> connectedComponenets) {

        resetBestConductance();

        for (int t : connectedComponenets.keySet()) {
            for (Set<Edge> edges : connectedComponenets.get(t)) {
                // For each set of edges in the connected components, use it as root node for a DAG
                constructDAG(edges, t, connectedComponenets);
            }
        }

        System.out.println("Printing best values");
        System.out.println("Conductance: "+bestConductance);
        System.out.println("Nodes: "+bestNodes);
        System.out.println("Start: "+bestStartTime+" End: "+bestEndTime);
        System.out.println("Zero Count: "+zeroCount);

        System.out.println("=============================================================");
    }

    public void constructDAG(Set<Edge> rootEdge, int startTime, TreeMap<Integer, List<Set<Edge>>> connectedComponenets ) {

        Set<Integer> rootNode = new HashSet<>(rootEdge.stream()
                .map(edge -> edge.toSet())
                .collect(HashSet::new, Set::addAll, Set::addAll));

        ArrayList<Set<Integer>> intersections = new ArrayList<>();
        intersections.add(rootNode);
        int currentTime = startTime;

        while(currentTime < (connectedComponenets.keySet().size() - 1)) {
            currentTime++;
            intersections = getIntersections(intersections, connectedComponenets.get(currentTime));
            if (intersections.size() != 0) {
                final int endTime = currentTime;
                intersections.parallelStream().forEach(nodes -> calculateConductance(nodes, startTime, endTime));
            } else {
                break;
            }

        }
    }

    public void calculateConductance(Set<Integer> nodes, int startTime, int endTime) {

//        System.out.println("Calculating conductance for size "+nodes.size()+" between "+startTime+" and "+endTime);

        int totalInternalWeight = 0;
        int totalExternalWeight = 0;
        int totalCutWeight = 0;


        for (int i=startTime;i<endTime;i++) {
            for (Edge e : timeGraph.get(i).keySet()) {
                int weight = timeGraph.get(i).get(e);
                if (e.listHasEdge(nodes)) totalInternalWeight += weight;
                else if (!e.listHasPartialEdge(nodes)) totalExternalWeight += weight;
                else  totalCutWeight += weight;
            }
        }
        totalExternalWeight += totalCutWeight;
        totalInternalWeight += totalCutWeight;

        double conductance = totalCutWeight / Math.min(totalInternalWeight, totalExternalWeight);
        if(conductance == 0) {
            zeroCount++;
        } else {
            // Normalize by root of (endtime - starttime)
            conductance = conductance / Math.sqrt(endTime - startTime);

            if (bestConductance == 0 || conductance < bestConductance) {
                bestConductance = conductance;
                bestNodes = nodes;
                bestStartTime = startTime;
                bestEndTime = endTime;
            }
        }
    }

    public ArrayList<Set<Integer>> getIntersections(List<Set<Integer>> parentNodes, List<Set<Edge>> cComps) {
        ArrayList<Set<Integer>> result = new ArrayList<>();
        for (Set<Edge> set : cComps) {
            for (Set<Integer> parent : parentNodes) {
                Set<Integer> intersection = new HashSet<>(set.stream()
                        .map(edge -> edge.toSet())
                        .collect(HashSet::new, Set::addAll, Set::addAll));
                intersection.retainAll(parent);

                if (intersection.size() > 1) { // Require intersection of size > 1
                    result.add(intersection);
                }
            }
        }
        return result;
    }

    public TreeMap<Integer, List<Set<Edge>>> getConnectedComponents(double k) {
        Set<Edge> thresholdedEdges = getThresholdedEdges(k);
        System.out.println("For k: "+k+" thresholded edges: "+thresholdedEdges.size());
        TreeMap<Integer, List<Set<Edge>>> result = new TreeMap<>();
        for (int time : timeGraph.keySet()) {
//            System.out.println("======================================================");
//            System.out.println("For time: " + time);

            HashMap<Edge, Integer> edgeMap = timeGraph.get(time);

            List<Set<Edge>> timeList = new ArrayList<Set<Edge>>();
            Set<Edge> connectedSet;

            Set<Edge> nodesToSearch = new HashSet<Edge>();
            nodesToSearch.addAll(edgeMap.keySet().parallelStream()
                    .filter(edge -> thresholdedEdges.contains(edge))
                    .collect(Collectors.toSet()));

            while (!nodesToSearch.isEmpty()) {
                Stack<Edge> stack = new Stack<Edge>();
                stack.push(nodesToSearch.stream().findFirst().get());
                connectedSet = new HashSet<Edge>();
                while (!stack.isEmpty()) {
                    Edge edge = stack.pop();
                    nodesToSearch.remove(edge);
                    Set<Edge> edgesToSearch = edgeMap.keySet().parallelStream()
                            .filter(e -> thresholdedEdges.contains(e)).collect(Collectors.toSet());

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
//            timeList.stream().forEach(set -> System.out.println("Set: " + set));
            result.put(time, timeList);
        }
        return result;
    }



    private void buildGraph(int node, int startTime, TreeMap<Integer, List<Set<Integer>>> input) {

        TreeMap<Integer, Set<Edge>> result = new TreeMap<>();
        Set<Integer> subSet = input.get(startTime).stream().filter(set -> set.contains(node)).findFirst().get();

        int time = startTime;
        while (time <= time + MAX_SEARCH_TIME) {
            Optional<Set<Integer>> setAtTime = input.get(time).stream().filter(set -> set.contains(node)).findFirst();
            if (setAtTime.isPresent()) {
                subSet.retainAll(setAtTime.get());
                time++;
            } else {
                break;
            }
        }
        if ((time - startTime) > 3) getConductance(result);
        else System.out.println("Too few components to be viable");
    }

    private void getConductance(TreeMap<Integer, Set<Edge>> dag) {
        int totalInternalWeights = 0;
        int totalExteralWeights = 0;
        int totalCutWeights = 0;
        for (int i : dag.keySet()) {
            Set<Edge> internalEdges = dag.get(i);
            Set<Edge> externalEdges = new HashSet<>();
            externalEdges.addAll(timeGraph.get(i).keySet());
            externalEdges.removeIf(edge -> edge.containsAny(internalEdges));
            Set<Edge> cutEdges = new HashSet<>();
            cutEdges.addAll(timeGraph.get(i).keySet());
            cutEdges.removeAll(internalEdges);
            cutEdges.retainAll(externalEdges);
            totalInternalWeights += timeGraph.get(i).keySet().stream()
                    .filter(e -> dag.get(i).contains(e))
                    .mapToInt(e -> timeGraph.get(i).get(e))
                    .sum();

            totalExteralWeights += timeGraph.get(i).keySet().stream()
                    .filter(edge -> externalEdges.contains(edge))
                    .mapToInt(value -> timeGraph.get(i).get(value))
                    .sum();
            totalCutWeights += timeGraph.get(i).keySet().stream()
                    .filter(edge -> cutEdges.contains(edge))
                    .mapToInt(value -> timeGraph.get(i).get(value))
                    .sum();
        }

        double conductance = totalCutWeights / Math.min(totalInternalWeights, totalExteralWeights);


    }

    private Set<Integer> getNeighbours(int time, int node) {
        Set<Integer> result = new HashSet<Integer>();
        if (timeGraph.containsKey(time)) {
            result.addAll(timeGraph.get(time).keySet().stream()
                    .filter(o -> o.containsNode(node))
                    .map(o -> o.getOtherNode(node)).collect(Collectors.toList()));
        }
        return result;
    }

    private void printGraph() {
        for (Integer i : timeGraph.keySet()) {
            System.out.println("Total Values for time: " + i + " is " + timeGraph.get(i).keySet());
        }
    }
}
