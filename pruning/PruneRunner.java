package pruning;

import java.io.File;
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
                    writer.println(l.size() + "," + l.stream().mapToInt(set -> set.size()).average().getAsDouble());
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
        while (((baseResult * 0.10) + baseResult) <= incResult) {
            incResult = getResult(prune, k + K_INCREMENT).normalizedConductance;

        }

        return baseResult;
    }

    public void printMatrix(String filename, double[][] matrix) {
        try {
            PrintWriter writer = new PrintWriter(filename);
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++) {
                    writer.write(matrix[i][j] + " ");
                }
                writer.write("\n");
            }
            writer.close();
        } catch (FileNotFoundException f) {
            System.err.println("Write Failed");
        }
    }

    public double getMetaInfo(String basepath, String filename, Prune prune) {
        double burstCond = 0;

        try {
            Scanner metaScanner = new Scanner(new File(basepath + "/meta/" + filename));

            // skip the first 6 lines to get to burst info
            int skip = 6;
            while (skip > 0) {
                metaScanner.nextLine();
                skip--;
            }
            Double startValue = Double.parseDouble(metaScanner.nextLine());
            Double endValue = Double.parseDouble(metaScanner.nextLine());
            int burstStart = startValue.intValue();
            int burstEnd = endValue.intValue();
            int burstCount = Integer.parseInt(metaScanner.nextLine());
            Set<Integer> bursNodes = new HashSet<>();
            for (int i = 0; i < burstCount; i++) {
                bursNodes.add(Integer.parseInt(metaScanner.nextLine()));
            }


            burstCond = prune.calculateConductance(bursNodes, burstStart, burstEnd);
            System.out.println("Burst between " + burstStart + " and " + burstEnd + " with " + burstCount + " nodes.");
            System.out.println("Raw Conductance of burst: " + burstCond);
            System.out.println("Normalized Conductance of burst: " +
                prune.normalizeByTime2(burstCond, burstStart, burstEnd));
            burstCond = prune.normalizeByTime2(burstCond, burstStart, burstEnd);
        } catch (FileNotFoundException fe) {
            fe.printStackTrace();
        }

        return burstCond;
    }

    private ArrayList<String> addIGapFiles(ArrayList<String> filenames, int range) {
        ArrayList<String> iGapFiles = new ArrayList<>();
        String fileStart = "iGap_series";

        iGapFiles.add(fileStart + "_020_run0");
        iGapFiles.add(fileStart + "_015_run0");
        iGapFiles.add(fileStart + "_010_run0");
        iGapFiles.add(fileStart + "_005_run0");
        iGapFiles.add(fileStart + "_002_run0");

        for (String name : iGapFiles) {
            for (int i = 0; i < range; i++) {
                filenames.add(name + i + ".txt");
            }
        }

        return filenames;
    }

    private ArrayList<String> addMiscFiles(ArrayList<String> filenames) {
        filenames.add("UBD080UBC040Gv001000Ge20Gg0100Bv0010Bd0005Bg01000400_10.txt");
        return filenames;
    }

    public static void main(String[] qwerty) {

        final int RANGE = 9;

        PruneRunner runner = new PruneRunner();

        ArrayList<String> filenames = new ArrayList<>();

//        filenames = runner.addMiscFiles(filenames);
//        filenames = runner.addIGapFiles(filenames, RANGE);
        filenames.add("iGap_series_002_run03.txt");

        String basepath = "output";
        double average_thresh_conductance = 0d;
        double average_burst_conductance = 0d;
        for (String filename : filenames) {

            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            System.out.println("Loading file:" + filename);
            Prune prune = runner.loadFile(filename, basepath);

            double burstCond = runner.getMetaInfo(basepath, filename, prune);
            average_burst_conductance += burstCond;

            /** Finding best conductance for k values **/
            double[] kList = {6};
//            for (double k : kList) {
//                System.out.println("For k: " + k);
//                TreeMap<Integer, List<Set<Edge>>> connectedComps = runner.getConnectedData(prune, k);
//                average_thresh_conductance += runner.getResult(prune, connectedComps).normalizedConductance;
//                prune.printTop5();
//                System.out.println("=============================================================");
//            }

            int prunedBounds;
            int n;
            for (int range = 0; range < 99; range++) {
                prunedBounds = 0;
                for (n = 0; n < (100 - range); n++) {
                    int end = n + range;
                    double result = prune.getBounds(n, end);
                    System.out.println(range+" "+n+" "+result);
                    if (result > burstCond) prunedBounds++;
                }
//                System.out.println(range+" "+n+" "+prunedBounds);
            }
            Instant start = Instant.now();
            Instant end = Instant.now();
            System.out.println("Duration: " + Duration.between(start, end).getSeconds());
        }

//        System.out.println(average_burst_conductance/10 + " "+average_thresh_conductance/10);
    }
}
