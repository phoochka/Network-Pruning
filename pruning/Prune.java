package pruning;


import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.IntToDoubleFunction;
import java.util.stream.Collectors;

//import org.la4j.decomposition.EigenDecompositor;
//import org.la4j.iterator.MatrixIterator;
//import org.la4j.matrix.SparseMatrix;
//import org.la4j.matrix.sparse.CRSMatrix;

public class Prune {


    private TreeMap<Integer, HashMap<Edge, Double>> timeGraph;
    //    private int totalNodes;
    private int maxNode;
    public int totalTime;

    public int startTime;
    public int endTime;

    private static final double epsilon = 0.1d;
    private static final int MIN_INTERSECTIONS = 1;

    private IntToDoubleFunction normalizingFunciton;

    private double bestNormConductance;
    private double bestRawConductance;
    private double bestN2factor;
    private Set<Integer> bestNodes;
    private int bestStartTime;
    private int bestEndTime;

    private SortedMap<Double, Result> results;


    public Prune(String filename, Path basepath, boolean isOld){
        if(isOld) readOldGraph(filename, basepath);
        else readGraph(filename, basepath);
        //printGraph();
        //System.out.println("All nodes at time 82: " + timeGraph.get(82).keySet());
        //System.out.println("All neighbours of 21 at time 82: " + getNeighbours(82, 21));
    }

    public void resetBestConductance() {
        bestNormConductance = bestRawConductance = 1.0d;
    }

    private void readOldGraph(String filename, Path basepath) {
        maxNode = Integer.MIN_VALUE;
        timeGraph = new TreeMap<Integer, HashMap<Edge, Double>>();
        Set<Integer> allNodes = new HashSet<>();
        // squashedGraph = new HashMap<Edge, Integer>();
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

                maxNode = node1 > maxNode ? node1 : maxNode;
                maxNode = node2 > maxNode ? node2 : maxNode;

                allNodes.add(node1);
                allNodes.add(node2);
                Double t = Double.parseDouble(s[2]);
                int floorInt = new Double(Math.floor(t)).intValue();
                Edge nodePair = new Edge(node1, node2);

                // squashedGraph.compute(nodePair, (k,v) -> v ==null ? 1 : v + 1);

                HashMap<Edge, Double> edgeMap;
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
                    edgeMap = new HashMap<Edge, Double>();
                    edgeMap.put(nodePair, 1.0d);
                    timeGraph.put(floorInt, edgeMap);
                }
            }

//            totalNodes = allNodes.stream().max(Comparator.naturalOrder()).get();
            startTime = timeGraph.firstKey();
            endTime = timeGraph.lastKey();
            totalTime = endTime - startTime;

            System.out.println("Max Nodes: " + maxNode);
            System.out.println("Total Time Length: " + totalTime);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // getStats();
    }

    private void readGraph(String filename, Path basepath) {
        maxNode = Integer.MIN_VALUE;
        timeGraph = new TreeMap<Integer, HashMap<Edge, Double>>();
        String line;

        try {
            Path filepath = basepath.resolve("new_graph");
            filepath = filepath.resolve(filename);
            Scanner timeGraphScanner = new Scanner(filepath.toFile());
            timeGraphScanner.nextLine(); // Skip heading line

            while (timeGraphScanner.hasNextLine()) {
                line = timeGraphScanner.nextLine();
                String s[] = line.split(",");
                if (s.length != 4) {
                    System.out.println("Input file in not in expected format");
                    System.out.println("Columns found: "+s.length+" in line: "+line);
                    timeGraphScanner.close();
                    throw new NoSuchElementException();
                }
                Integer node1 = Integer.parseInt(s[0]);
                Integer node2 = Integer.parseInt(s[1]);
                Integer time = Integer.parseInt(s[2]);
                Double weight = Double.parseDouble(s[3]);

                maxNode = node1 > maxNode ? node1 : maxNode;
                maxNode = node2 > maxNode ? node2 : maxNode;

                Edge edge = new Edge(node1, node2);

                if (timeGraph.containsKey(time)) {
                    timeGraph.get(time).put(edge, weight);
                } else {
                    HashMap<Edge, Double> edgeMap = new HashMap<Edge, Double>();
                    edgeMap.put(edge, weight);
                    timeGraph.put(time, edgeMap);
                }
            }

            startTime = timeGraph.firstKey();
            endTime = timeGraph.lastKey();
            totalTime = endTime - startTime;

//            System.out.println("Max Nodes: " + maxNode);
            System.out.println("Total Time Length: " + totalTime);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public TreeMap<Integer, HashMap<Integer, Double>> getNodeWeights() {
        TreeMap<Integer, HashMap<Integer, Double>> nodeGraph = new TreeMap<>();

        for (int time : timeGraph.keySet()) {
            if (!nodeGraph.containsKey(time)) nodeGraph.put(time, new HashMap<Integer, Double>());
            for (Map.Entry<Edge, Double> entry : timeGraph.get(time).entrySet()) {
                nodeGraph.get(time).compute(entry.getKey().getNode1(), (k, v) -> v == null ? entry.getValue() : v + entry.getValue());
                nodeGraph.get(time).compute(entry.getKey().getNode2(), (k, v) -> v == null ? entry.getValue() : v + entry.getValue());
            }
        }

        return nodeGraph;
    }

    public double getIntervalWeight(TreeMap<Integer, HashMap<Integer, Double>> nodeGraph, int subStart, int subEnd, int wholeStart, int wholeEnd) {

        double weight = Double.MAX_VALUE;
        for (int node : getNodesInTime(subStart, subEnd)) {
            double ratio = getDegreeForNode(nodeGraph, subStart, subEnd, node) / getDegreeForNode(nodeGraph, wholeStart, wholeEnd, node);
            weight = ratio < weight ? ratio : weight;
        }

        return weight;
    }

    private double getDegreeForNode(TreeMap<Integer, HashMap<Integer, Double>> nodeGraph, int start, int end, int node) {
        return nodeGraph.subMap(start, true, end, true).values().stream()
            .filter(keys -> keys.containsKey(node))
            .mapToDouble(map -> map.get(node)).sum();
    }

    private double getDegreeForNode(int start, int end, int node) {
        double sum = 0;
        for (int i = start; i <= end; i++) {
            sum += timeGraph.get(i).entrySet().stream()
                .filter(e -> e.getKey().containsNode(node))
                .mapToDouble(e -> e.getValue())
                .sum();
        }
        return sum;
    }


    public Set<Integer> getNodesInTime(int start, int end) {
        Set<Integer> nodes = new HashSet<>();
        for (int i = start; i <= end; i++) {
            timeGraph.get(i).keySet().forEach(k -> nodes.addAll(k.toSet()));
        }

        return nodes;
    }


    public double getTotalWeights(int start, int end) {
        double sum = 0;
        for (int i = start; i <= end; i++) {
            sum += timeGraph.get(i).values().stream().mapToDouble(Double::doubleValue).sum();
        }
        return sum;
    }


    public double getBounds(int start, int end) {

        LinkedSparseMatrix lsm = new LinkedSparseMatrix(maxNode+1, maxNode+1);

        for (int i = start; i <= end; i++) {
            for (Edge e : timeGraph.get(i).keySet()) {
                int n1 = e.getNode1();
                int n2 = e.getNode2();
                double weight = timeGraph.get(i).get(e);
                lsm.set(n1, n2, lsm.get(n1,n2) - weight);
                lsm.set(n2, n1, lsm.get(n2,n1) - weight);

                lsm.set(n1, n1, lsm.get(n1,n1) + weight);
                lsm.set(n2, n2, lsm.get(n2,n2) + weight);
            }
        }

        LinkedSparseMatrix normalizedLaplacian = new LinkedSparseMatrix(maxNode+1, maxNode+1);

        Iterator<MatrixEntry> iter = lsm.iterator();
        while(iter.hasNext()) {
            MatrixEntry me = iter.next();
            int c = me.column();
            int r = me.row();
            double val = me.get();

            double normalizedVal = val / (Math.sqrt(lsm.get(c,c)) * Math.sqrt(lsm.get(r,r)));

            normalizedLaplacian.set(r,c,normalizedVal);
        }

        ArpackSym solver = new ArpackSym(normalizedLaplacian);
        Map<Double, DenseVectorSub> resultsSA = solver.solve((int)Math.log(maxNode), ArpackSym.Ritz.SM);

        double lambda2 = (double) resultsSA.keySet().toArray()[resultsSA.keySet().size()-2];
        double lowerBounds = lambda2 / 2.0;

        return lowerBounds;

    }

    public Set<Edge> thresholdEdgesByPer(double k, HashMap<Edge, Double> edgeWeights) {
        int edgeListSize = (int) (edgeWeights.size() * (k / 100));
        Comparator<Map.Entry<Edge, Double>> sortEdgesByWeight = (e1, e2) -> e1.getValue().compareTo(e2.getValue());
        return edgeWeights.entrySet().stream()
            .sorted(sortEdgesByWeight.reversed())
            .limit(edgeListSize)
            .map(e -> e.getKey())
            .collect(Collectors.toSet());
    }

    public Result traverseGraph(TreeMap<Integer, List<Set<Edge>>> connectedComponents, IntToDoubleFunction f) {

        this.normalizingFunciton = f;

        resetBestConductance();

        int dagCount = 0;

        //results = Collections.synchronizedSortedMap(new TreeMap<Double, Result>());
        results = new ConcurrentSkipListMap<Double,Result>();

        for (int t : connectedComponents.keySet()) {
            for (Set<Edge> edges : connectedComponents.get(t)) {
                // For each set of edges in the connected components, use it as root node for a DAG
                constructDAG(edges, t, connectedComponents);
                dagCount++;
            }
        }
//        System.out.println("DAG COUNT: " + dagCount);
        Result result = new Result(bestNormConductance, bestRawConductance, bestNodes, bestStartTime, bestEndTime);

        System.out.println("Normalized Conductance: " + bestNormConductance);
        System.out.println("Raw ConductanceL " + bestRawConductance);
        System.out.println("N2 Factor: " + bestN2factor);
//        System.out.println("Nodes: "+bestNodes);
        System.out.println("No. of nodes: " + bestNodes.size());
        System.out.println("Start: " + bestStartTime + " End: " + bestEndTime);

//        printTraversalStats(traversals);

        return result;
    }

    private void printTraversalStats(List<Integer> traversals) {
        traversals.stream().forEach(i -> System.out.print(i + ", "));
        System.out.println("\n");
        HashMap<Integer, Integer> traversalCounter = new HashMap<>();
        traversals.stream().forEach(i -> traversalCounter.merge(i, 1, Math::addExact));
        traversalCounter.entrySet().stream()
            .sorted(Map.Entry.<Integer, Integer>comparingByValue().reversed())
            .forEach(k -> System.out.print(k + " "));

    }

    public void constructDAG(Set<Edge> rootEdge, int startTime, TreeMap<Integer, List<Set<Edge>>> connectedComponenets) {

        Set<Integer> rootNode = new HashSet<Integer>();
        for (Edge e : rootEdge) {
            rootNode.addAll(e.toSet());
        }
        ArrayList<Set<Integer>> intersections = new ArrayList<>();
        intersections.add(rootNode);
        int currentTime = startTime;
        int traversalCount = 0;
        while (currentTime < connectedComponenets.lastKey()) {
            traversalCount++;

            if (connectedComponenets.get(++currentTime) == null) break;
            else {
                intersections = getIntersections(intersections, connectedComponenets.get(currentTime));
                if (intersections.size() != 0) {
                    final int endTime = currentTime;
                    intersections.parallelStream().forEach(nodes -> updateConductance(nodes, startTime, endTime));
                } else break;
            }
        }
    }

    public void updateConductance(Set<Integer> nodes, int startTime, int endTime) {

        double rawConductance = calculateConductance(nodes, startTime, endTime);

        if (rawConductance < 0.000001d) {
            // TODO: HANDLE ZERO CONDUCTANCE
            System.err.println("Raw cond is 0");
        } else {
            // Normalize
            double normalizedConductance = normalize(rawConductance, startTime, endTime);

            // Store Top K
            Result r = new Result(normalizedConductance, rawConductance, nodes, startTime, endTime);
            updateResults(rawConductance, r);

            if (normalizedConductance < bestNormConductance) {
//                if (rawConductance < bestRawConductance) {

                bestNormConductance = normalizedConductance;
                bestRawConductance = rawConductance;
                bestN2factor = getN2factor(startTime, endTime);
                bestNodes = nodes;
                bestStartTime = startTime;
                bestEndTime = endTime;
            }
        }
    }



    private void updateResults(Double conductance, Result r) {

        // Checking if this time period already exists
//        final SortedMap<Double, Result> checkMap = new TreeMap<Double, Result>();
        if (results.size() != 0 ) {
            for (Result storedResult : results.values()) {
                if (conductance > r.normalizedConductance && r.startTime > storedResult.startTime && r.endTime < storedResult.endTime)
                    return;
            }
        }

        results.put(r.normalizedConductance, r);

        if (results.size() > 100) {
            Double lastKey = results.lastKey();
            results.remove(lastKey);
        }
    }

    public double calculateConductance(Set<Integer> nodes, int startTime, int endTime) {

//        System.out.println("Calculating conductance for size "+nodes.size()+" between "+startTime+" and "+endTime);

        double totalInternalWeight = 0d;
        double totalExternalWeight = 0d;
        double totalCutWeight = 0d;

        for (int i = startTime; i < endTime; i++) {
            for (Edge e : timeGraph.get(i).keySet()) {
                double weight = timeGraph.get(i).get(e);
                if (e.listHasEdge(nodes)) totalInternalWeight += weight;
                else if (e.listHasPartialEdge(nodes)) totalCutWeight += weight;
                else totalExternalWeight += weight;
            }
        }
        totalExternalWeight = totalExternalWeight + totalCutWeight;
        totalInternalWeight = totalInternalWeight + totalCutWeight;

        double conductance = totalCutWeight / Math.min(totalInternalWeight, totalExternalWeight);
        return conductance;
    }

    private double normalize(double rawConductance, int startTime, int endTime) {
        double normalizeBy = normalizingFunciton.applyAsDouble(endTime - startTime);
        return rawConductance*normalizeBy;
    }

    public double normalizeByTime1(double conductance, int startTime, int endTime) {
        return conductance / Math.sqrt(endTime - startTime);
    }

    public double normalizeByTime2(double conductance, int startTime, int endTime) {
        return conductance / Math.cbrt(endTime - startTime);
    }

    public double getN2factor(int startTime, int endTime) {
        return Math.cbrt(endTime - startTime);
    }

    public double normalizeByLog(double conductance, int startTime, int endTime) {
        return conductance / Math.log(endTime - startTime);
    }

    public double normalizeByE(double e, double conductance, int startTime, int endTime) {
        double d = endTime - startTime;
        return conductance / Math.exp(d * e);
    }

    public ArrayList<Set<Integer>> getIntersections(List<Set<Integer>> parentNodes, List<Set<Edge>> cComps) {
        ArrayList<Set<Integer>> result = new ArrayList<>();

        for (Set<Edge> set : cComps) {
            for (Set<Integer> parent : parentNodes) {
//                Set<Integer> intersection = new HashSet<>(set.stream()
//                        .map(edge -> edge.toSet())
//                        .collect(HashSet::new, Set::addAll, Set::addAll));

                Set<Integer> intersection = new HashSet<Integer>();
                for (Edge e : set) {
                    intersection.addAll(e.toSet());
                }

                intersection.retainAll(parent);

                if (intersection.size() >= MIN_INTERSECTIONS) { // Require intersection of size > 1
                    result.add(intersection);
                }
            }
        }
        return result;
    }

    public void checkConnectedComponenets(int time) {
        HashMap<Edge, Double> edgeMap = timeGraph.get(time);
        List<Set<Edge>> timeList = new ArrayList<Set<Edge>>();
        Set<Edge> connectedSet;

        Set<Edge> nodesToSearch = new HashSet<Edge>();
        nodesToSearch.addAll(edgeMap.keySet().parallelStream()
            .collect(Collectors.toSet()));

        while (!nodesToSearch.isEmpty()) {
            Stack<Edge> stack = new Stack<Edge>();
            stack.push(nodesToSearch.stream().findFirst().get());
            connectedSet = new HashSet<Edge>();
            while (!stack.isEmpty()) {
                Edge edge = stack.pop();
                nodesToSearch.remove(edge);
                Set<Edge> edgesToSearch = edgeMap.keySet().parallelStream()
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

        System.out.println("Number of connected comps at time: " + time + " = " + timeList.size());
        timeList.stream().filter(s -> s.size() < 10).forEach(s -> System.out.println(s));
    }

    public TreeMap<Integer, List<Set<Edge>>> getConnectedComponents(double k) {

        TreeMap<Integer, List<Set<Edge>>> result = new TreeMap<>();
        for (int time : timeGraph.keySet()) {

            HashMap<Edge, Double> edgeMap = timeGraph.get(time);

            Set<Edge> thresholdedEdges = thresholdEdgesByPer(k, edgeMap);

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
                        .filter(thresholdedEdges::contains).collect(Collectors.toSet());

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
            result.put(time, timeList);
        }
        return result;
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

    public void printGraph() {
        for (Integer i : timeGraph.keySet()) {
            System.out.println("Total Values for time: " + i + " is " + timeGraph.get(i).keySet());
        }
    }


    public void printTop(int k) {
        System.out.println("PRINTING TOP "+k);
        if (results == null) System.out.println("NO RESULTS");
        results.keySet().stream().limit(k)
            .forEach(c -> System.out.println(
                      "Raw Cond: "   + results.get(c).rawConductance
                    + " Norm Cond: " + results.get(c).normalizedConductance
                    + " Size: "      + results.get(c).nodes.size()
                    + " Time: "      + results.get(c).startTime + " " + results.get(c).endTime
            ));
    }

}
