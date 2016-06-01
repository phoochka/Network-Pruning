package lkml_reply;


import pruning.Edge;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by gaurav on 5/30/16.
 */
public class ProcessMailingListData {

    final static long EPOCH_MONTH = 2629743;


    public ProcessMailingListData () {
        TreeMap<Integer, HashMap<Edge, Integer>> timeGraph = processFile(new File("lkml_reply/out.lkml-reply"));
        connectingGraph(timeGraph);
        writeToFile(timeGraph);
    }

    private void writeToFile(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {
        try {
            PrintWriter writer = new PrintWriter("output/new_graph/mailing_list_data.txt");
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

        for (int time : timeGraph.keySet()){
            for (Edge edge : allEdges){
                if (!timeGraph.get(time).keySet().contains(edge)) timeGraph.get(time).put(edge, 1);
            }
        }
        System.out.println("Finished connecting graph");
    }


    private TreeMap<Integer, HashMap<Edge, Integer>> processFile(File f) {

        TreeMap<Integer, HashMap<Edge, Integer>> timeGraph = new TreeMap<>();

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

                    long epochTime = Long.parseLong(line[3]);
                    int time = Math.toIntExact(epochTime / (EPOCH_MONTH*24));


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

        System.out.println("Finished creating data structure");

        return timeGraph;
    }

    public static void main (String[] sdff) {
        new ProcessMailingListData();
    }
}

