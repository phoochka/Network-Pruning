package pruning;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

import org.jblas.*;


class Prune {


    private TreeMap<Integer, HashMap<Edge, Double>> timeGraph;
    private int totalNodes;

    private static final double epsilon = 0.1d;
    private static final int MIN_INTERSECTIONS = 2;

    // private HashMap<Edge, Integer> squashedGraph;
    // private int size;
    // private double mean;
    // private double stdDev;

    private double bestNormConductance;
    private double bestRawConductance;
    private double bestN2factor;
    private Set<Integer> bestNodes;
    private int bestStartTime;
    private int bestEndTime;

    /****
     * TESTING
     *****/
    public int zeroConductanceCount;
    public int calculateCondCount;
    public List<Integer> traversals;
    /******************/

    private SortedMap<Double, Result> results;


    public void resetBestConductance() {
        bestNormConductance = 1.0d;
        zeroConductanceCount = 0;
        calculateCondCount = 0;
        traversals = Collections.synchronizedList(new ArrayList<>());
    }

    public Prune(String filename, String basepath) {
        readGraph(filename, basepath);
        //printGraph();
        //System.out.println("All nodes at time 82: " + timeGraph.get(82).keySet());
        //System.out.println("All neighbours of 21 at time 82: " + getNeighbours(82, 21));
    }

    private void readGraph(String filename, String basepath) {
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

            totalNodes = allNodes.size();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // getStats();
    }

    public double[][] matrixMultiply(double[][] A, double[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (nA != mB) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] C = new double[mA][nB];
        for (int i = 0; i < mA; i++)
            for (int j = 0; j < nB; j++)
                for (int k = 0; k < nA; k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    public double[][] testMatrix() {
        double[][] degreeMatrix = {{3, 0, 0}, {0, 4, 0}, {0, 0, 5}};
        double[][] laplacian = {{3, -1, -2}, {-1, 4, -3}, {-2, -3, 5}};

        // Getting inverse square root of the diagonal
        for (int n = 0; n < degreeMatrix.length; n++) {
            System.out.println("Changing value of " + degreeMatrix[n][n]);
            degreeMatrix[n][n] = 1 / (Math.sqrt(degreeMatrix[n][n]));
        }

        for (int i = 0; i < degreeMatrix.length; i++) {
            for (int j = 0; j < degreeMatrix.length; j++) {
                System.out.print(degreeMatrix[i][j] + " ");
            }
            System.out.println();
        }

        double[][] result = matrixMultiply(matrixMultiply(degreeMatrix, laplacian), degreeMatrix);

        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result.length; j++) {
                System.out.print(result[i][j] + " ");
            }
            System.out.println();
        }
        return result;
    }

    public double[][] getNormalizedLaplacianStupidly(int start, int end) {

        double[][] laplacian = new double[totalNodes][totalNodes];
        double[][] degreeMatrix = new double[totalNodes][totalNodes];
        for (int i = start; i < end; i++) {
            for (Edge e : timeGraph.get(i).keySet()) {

                double weight = timeGraph.get(i).get(e);

                laplacian[e.getNode1()][e.getNode2()] -= weight;
                laplacian[e.getNode2()][e.getNode1()] -= weight;


                laplacian[e.getNode1()][e.getNode1()] += weight;
                laplacian[e.getNode2()][e.getNode2()] += weight;

                degreeMatrix[e.getNode1()][e.getNode1()] += weight;
                degreeMatrix[e.getNode2()][e.getNode2()] += weight;
            }
        }

        // Getting inverse square root of the diagonal
        for (int n = 0; n < totalNodes; n++) {
            degreeMatrix[n][n] = 1 / Math.sqrt(degreeMatrix[n][n]);
        }

        double[][] normalizedLaplacian = matrixMultiply(matrixMultiply(degreeMatrix, laplacian), degreeMatrix);
        return normalizedLaplacian;
    }

    public double[][] getNormalizedLaplacian(int start, int end) {
        double[][] laplacian = new double[totalNodes][totalNodes];
        for (int i = start; i <= end; i++) {
            for (Edge e : timeGraph.get(i).keySet()) {
                double weight = timeGraph.get(i).get(e);

                laplacian[e.getNode1()][e.getNode2()] -= weight;
                laplacian[e.getNode2()][e.getNode1()] -= weight;

                laplacian[e.getNode1()][e.getNode1()] += weight;
                laplacian[e.getNode2()][e.getNode2()] += weight;
            }
        }

        /** ADD A DEFAULT WEIGHT OF 0.1 **/
        for (int i = 0; i < totalNodes; i++) {
                laplacian[i][i] += 0.1;
        }

        double[][] normalizedLaplacian = new double[totalNodes][totalNodes];
        for (int i = 0; i < totalNodes; i++) {
            for (int j = 0; j < totalNodes; j++) {

                if (laplacian[i][i] == 0.0) laplacian[i][i] = 0.1d;
                if (laplacian[j][j] == 0.0) laplacian[j][j] = 0.1d;

                normalizedLaplacian[i][j] = laplacian[i][j] / (Math.sqrt(laplacian[i][i]) * Math.sqrt(laplacian[j][j]));
                if (i == j && Double.isNaN(normalizedLaplacian[i][j])) {
                    System.out.println("For " + start + " " + end + " At " + i + "," + j + " NAN: " + laplacian[i][j] + " " + laplacian[i][i] + " " + laplacian[j][j]);
                    normalizedLaplacian[i][j] = 1.0d;
                } else if (i == j && laplacian[i][j] == 0.0d) System.out.println("WTF");
            }
        }

        return normalizedLaplacian;
    }

    public double getBounds(int start, int end) {
        double result = 0;

        DoubleMatrix A = new DoubleMatrix(getNormalizedLaplacian(start, end));
        DoubleMatrix B = DoubleMatrix.eye(A.rows);
        try {
            DoubleMatrix La = Eigen.symmetricGeneralizedEigenvalues(A, B, 1, 1);
            double lambda2 = La.get(0, 0);
            result = lambda2 / 2;
        } catch (org.jblas.exceptions.LapackConvergenceException lce) {
            System.err.println("Convergence Error");
            result = -1;
        } catch (Exception e) {
            System.err.println("Error "+e);
        }
        return result;
    }

    /**
     * private void getStats() {
     * <p>
     * size = squashedGraph.values().size();
     * <p>
     * // Finding mean
     * // int sum = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).sum();
     * mean = squashedGraph.values().stream().mapToInt(i -> new Integer(i)).average().getAsDouble();
     * <p>
     * // System.out.println("Mean is: "+mean);
     * <p>
     * // Finding standard deviation
     * double temp = 0;
     * for (int i : squashedGraph.values()) {
     * temp += (mean - i) * (mean - i);
     * }
     * stdDev = Math.sqrt(temp / size);
     * // System.out.println("stddev is: "+stdDev);
     * }
     * <p>
     * <p>
     * public Set<Edge> thresholdEdgesByStd(double k) {
     * double thresh = mean + (k * stdDev);
     * return squashedGraph.keySet().parallelStream()
     * .filter(e -> squashedGraph.get(e) > thresh).collect(Collectors.toSet());
     * }
     **/

    public Set<Edge> thresholdEdgesByPer(double k, HashMap<Edge, Double> edgeWeights) {
        int edgeListSize = (int) (edgeWeights.size() * (k / 100));
        Comparator<Map.Entry<Edge, Double>> sortEdgesByWeight = (e1, e2) -> e1.getValue().compareTo(e2.getValue());
        return edgeWeights.entrySet().stream()
            .sorted(sortEdgesByWeight.reversed())
            .limit(edgeListSize)
            .map(e -> e.getKey())
            .collect(Collectors.toSet());
    }

    public Result traverseGraph(TreeMap<Integer, List<Set<Edge>>> connectedComponents) {

        resetBestConductance();

        int dagCount = 0;

        results = Collections.synchronizedSortedMap(new TreeMap<Double, Result>());

        for (int t : connectedComponents.keySet()) {
            for (Set<Edge> edges : connectedComponents.get(t)) {
                // For each set of edges in the connected components, use it as root node for a DAG
                constructDAG(edges, t, connectedComponents);
                dagCount++;
            }
        }
        System.out.println("DAG COUNT: " + dagCount);
        System.out.println("CALCULATE COND COUNT: " + calculateCondCount);
        Result result = new Result(bestNormConductance, bestRawConductance, bestNodes, bestStartTime, bestEndTime);

        System.out.println("Printing best values");
        System.out.println("Normalized Conductance: " + bestNormConductance);
        System.out.println("Raw ConductanceL "+bestRawConductance);
        System.out.println("N2 Factor: "+bestN2factor);
//        System.out.println("Nodes: "+bestNodes);
        System.out.println("No. of nodes: " + bestNodes.size());
        System.out.println("Start: " + bestStartTime + " End: " + bestEndTime);
        System.out.println("No. of zero cond results: " + zeroConductanceCount);

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
        traversals.add(traversalCount);
    }

    public void updateConductance(Set<Integer> nodes, int startTime, int endTime) {
        calculateCondCount++;

        double rawConductance = calculateConductance(nodes, startTime, endTime);

        if (rawConductance < 0.000001d) {
            zeroConductanceCount++;
        } else {
            // Normalize
            double normalizedConductance = normalizeByTime2(rawConductance, startTime, endTime);

            // Store Top K
            Result r = new Result(normalizedConductance, rawConductance, nodes, startTime, endTime);
            updateResults(normalizedConductance, r);

            if (normalizedConductance < bestNormConductance) {
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
        results.put(conductance, r);

        if (results.size() > 5) {
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

    private void printGraph() {
        for (Integer i : timeGraph.keySet()) {
            System.out.println("Total Values for time: " + i + " is " + timeGraph.get(i).keySet());
        }
    }


    public void printTop5() {
        System.out.println("PRINTING TOP 5");
        if (results == null) System.out.println("NO RESULTS");
        results.keySet().stream().limit(5)
            .forEach(k -> System.out.println(
                "Raw Cond: " + results.get(k).rawConductance
                    + "Norm Cond: " + results.get(k).normalizedConductance
                    + " Size: " + results.get(k).nodes.size()
                    + " Time: " + results.get(k).startTime + " " + results.get(k).endTime
            ));
    }

}
