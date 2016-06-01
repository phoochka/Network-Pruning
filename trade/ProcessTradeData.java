package trade;

import pruning.Edge;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by gaurav on 5/26/16.
 */
public class ProcessTradeData {

    TreeMap<Integer, HashMap<Edge, Integer>> timeGraph = new TreeMap<>();
    List<String> countries = new ArrayList<>();

    public ProcessTradeData() {
        File folder = new File("./trade/total_graphs");
        File[] listOfFiles = folder.listFiles();

        Arrays.stream(listOfFiles)
            .filter(thing -> thing.isFile())
            .filter(file -> file.getName().endsWith("_0.txt"))
            .forEach(f -> processFile(f));

        connectingGraph(timeGraph);
        writeToFile(timeGraph);
        writeToLogFile(timeGraph);
//        writeMapping(countries);
    }

    private void connectingGraph(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {

        Set<Edge> allEdges = new HashSet<>();
        timeGraph.values().parallelStream().map(map -> map.keySet()).forEach(e -> allEdges.addAll(e));

        for (int time : timeGraph.keySet()){
            for (Edge edge : allEdges){
                if (!timeGraph.get(time).keySet().contains(edge)) timeGraph.get(time).put(edge, 1);
            }
        }
    }

    private void writeMapping(List<String> countries) {
        try{
            PrintWriter writer = new PrintWriter("trade_data_mapping.txt");
            System.out.println("Writing countries of size "+countries.size());
            for (String s : countries) {
                writer.println(countries.indexOf(s)+" : "+s);
            }
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void writeToFile(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {
        try {
            PrintWriter writer = new PrintWriter("trade_data.txt");
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

    private void writeToLogFile(TreeMap<Integer, HashMap<Edge, Integer>> timeGraph) {
        try {
            PrintWriter writer = new PrintWriter("trade_data_log.txt");
            writer.println("Node1,Node2,Time,Weight");
            for (int time : timeGraph.keySet()) {
                for (Edge edge : timeGraph.get(time).keySet()) {
                    Integer weight = timeGraph.get(time).get(edge);
                    int int_weight = (int) (1.0d + Math.log(weight));
                    writer.println(edge.getNode1() + "," + edge.getNode2() + "," + time + "," + int_weight);
                }
            }
//            writer.println("FINISHED");
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    void processFile(File f) {

        int time = Integer.parseInt(f.getName().substring(0,4));
        String[] line;
        try {
            Scanner rawScanner = new Scanner(f);
            while(rawScanner.hasNextLine()){
               line = rawScanner.nextLine().split("\t");
                if (!countries.contains(line[0])) countries.add(line[0]);
                if (!countries.contains(line[1])) countries.add(line[1]);

                int node1 = countries.indexOf(line[0]);
                int node2 = countries.indexOf(line[1]);

                Edge edge = new Edge(node1, node2);

                int weight = Integer.parseInt(line[2]);

                if (timeGraph.containsKey(time)) {
                    if (timeGraph.get(time).containsKey(edge)) {
                        int oldWeight = timeGraph.get(time).get(edge);
                        timeGraph.get(time).replace(edge, oldWeight, weight+oldWeight);
                    } else timeGraph.get(time).put(new Edge(node1, node2), weight);
                } else {
                    HashMap<Edge, Integer> edgeWeights = new HashMap<>();
                    edgeWeights.put(edge, weight);
                    timeGraph.put(time, edgeWeights);
                }
            }
            rawScanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void main (String[] eroighe) {
        new ProcessTradeData();
    }
}
