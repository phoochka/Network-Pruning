package pruning;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.time.Duration;
import java.time.Instant;
import java.util.*;

public class PruneRunner {

    public PruneRunner() {
    }

    public Prune loadFile(String filename, String basepath) {
        return new Prune(filename, basepath);
    }

    public TreeMap<Integer, List<Set<Edge>>> getConnectedData(Prune prune, double k) {
        return prune.getConnectedComponents(k);
    }

    public Result getResult(Prune prune, TreeMap<Integer, List<Set<Edge>>> connectedComps) {
        return prune.traverseGraph(connectedComps);
    }

    public Result getResult(Prune prune, double k) {
        return prune.traverseGraph(prune.getConnectedComponents(k));
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
            System.err.println("Write Failed");
        }
    }

    public double thresholdBestResult(String filename, String basepath) {
        Prune prune = loadFile(filename, basepath);
        final double STARTING_K = 1.0d;
        final double ACCEPTABLE_CONDUCTANCE = 10;
        final double K_INCREMENT = 1;

        double baseResult = getResult(prune, STARTING_K).normalizedConductance;
        double incResult = Double.MAX_VALUE;
        double k = STARTING_K;
        while(((baseResult * 0.10)+baseResult) <= incResult) {
            incResult = getResult(prune, k+K_INCREMENT).normalizedConductance;

        }

        return baseResult;
    }

    private ArrayList<String> addIGapFiles(ArrayList<String> filenames) {
        String fileStart = "iGap";
        filenames.add(fileStart+"_series_002_run0");
        filenames.add(fileStart+"_series_005_run0");
        filenames.add(fileStart+"_series_020_run0");
        filenames.add(fileStart+"_series_015_run0");
        filenames.add(fileStart+"_series_010_run0");



        return filenames;
    }

    private ArrayList<String> addMiscFiles(ArrayList<String> filenames){
        filenames.add("UBD080UBC040Gv001000Ge20Gg0100Bv0010Bd0005Bg01000400_10");
        return  filenames;
    }

    public static void main(String[] qwerty) {

        PruneRunner runner = new PruneRunner();
        // String filename = "UBD080UBC040Gv001000Ge20Gg0100Bv0010Bd0005Bg01000400_10.txt";
        ArrayList<String> filenames = new ArrayList<>();

//        filenames = runner.addMiscFiles(filenames);
        filenames = runner.addIGapFiles(filenames);

//        filenames.add(fileStart+"_series_070_010_run0");
//        filenames.add(fileStart+"_series_070_015_run0");
        String basepath = "output";
        for (String filename : filenames) {
            double average_conductance = 0d;
            for (int i = 0; i < 1; i++) {
                String current_filename = filename+i+".txt";
                System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                System.out.println("Loading file:"+current_filename);
                Prune prune = runner.loadFile(current_filename, basepath);

                //outputKFactor(prune, kList);

                /** Finding synthetic burst data **/
                /**
                 Set<Integer> burstNodes = new HashSet<>(Arrays.asList(
                 107,110,134,149,216,251,443,477,901,996
                 ));
                 int burstStart = 20;
                 int burstEnd = 25;
                 double burstConductance = prune.calculateConductance(burstNodes, burstStart, burstEnd);
                 System.out.println("Burst Conductance: "+burstConductance);
                 System.out.println("Normalized1 conductance "+prune.normalizeByTime1(burstConductance, burstStart, burstEnd));
                 System.out.println("Normalized2 conductance "+prune.normalizeByTime2(burstConductance, burstStart, burstEnd));
                 System.out.println("Norm E 0.0175 conductance "+prune.normalizeByE(0.0175,burstConductance, burstStart, burstEnd));
                 System.out.println("Norm by Log conductance "+prune.normalizeByLog(burstConductance, burstStart, burstEnd));
                 **/

                /** Finding best conductance for k values **/
//                double[] kList = {1, 3, 6, 10};
                double[] kList = {20};
                for (double k: kList) {
                    System.out.println("For k: "+k);
                    for(int n=0;n<98;n+=2){
                        Instant start = Instant.now();
                        System.out.println(prune.getBounds(n, n+2)+"for "+n);
                        Instant end = Instant.now();
                        System.out.println("Duration: "+ Duration.between(start,end).getSeconds());
                    }
                    // TreeMap<Integer, List<Set<Edge>>> connectedComps = runner.getConnectedData(prune, k);
                    // runner.getResult(prune, connectedComps);
                    // prune.printTop5();
                    System.out.println("=============================================================");
                }
            }
            // Getting Average
        }

    }
}
