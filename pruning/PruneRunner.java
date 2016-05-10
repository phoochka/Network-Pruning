package pruning;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.time.Duration;
import java.time.Instant;
import java.util.*;

public class PruneRunner {

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


    public Result getMetaInfo(String basepath, String filename, Prune prune) {
        Result result = null;

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


            double rawBurstCond = prune.calculateConductance(bursNodes, burstStart, burstEnd);
            double normBurstCond = prune.normalizeByTime2(rawBurstCond, burstStart, burstEnd);
            System.out.println("Burst between " + burstStart + " and " + burstEnd + " with " + burstCount + " nodes.");
            System.out.println("Raw Conductance of burst: " + rawBurstCond);
            System.out.println("Normalized Conductance of burst: " + normBurstCond);

            result = new Result(rawBurstCond, normBurstCond, bursNodes, burstStart, burstEnd);
        } catch (FileNotFoundException fe) {
            fe.printStackTrace();
        }

        return result;
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


    public Result doThresholds(String basepath, String filename, Prune prune) {
        System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        System.out.println("Loading file:" + filename);
        Result result = null;
        double cond = -1;
        /** Finding best conductance for k values **/
        double[] kList = {1};
        for (double k : kList) {
            System.out.println("For k: " + k);
            TreeMap<Integer, List<Set<Edge>>> connectedComps = getConnectedData(prune, k);
            result = getResult(prune, connectedComps);
            prune.printTop5();
            System.out.println("=============================================================");
        }

        return result;
    }

    public void printBounds(Prune prune) {

        int n;
        int total = 0;
        for (int range = 0; range < 100; range++) {
            Instant tStart = Instant.now();
            for (n = 0; n < (100 - range); n++) {
//                if (range == 15 && n == 0) n = 14;
                int end = n + range;
                total++;
                System.out.println(n + " " + end + " " + prune.getBounds(n, end));
            }
        }
        System.out.println("Total: " + total);
    }

    public HashMap processesBounds1(String filename) {
        HashMap<Integer, TreeMap<Integer, Double>> allBounds = new HashMap<>();
        try {
            Scanner boundScanner = new Scanner(new File(filename));
            while (boundScanner.hasNextLine()) {
                String values[] = boundScanner.nextLine().split(" ");
                int timeLength = Integer.parseInt(values[0]);
                int index = Integer.parseInt(values[1]);
                double bound = Double.parseDouble(values[2]);
                if (allBounds.containsKey(timeLength)) {
                    allBounds.get(timeLength).put(index, bound);
                } else {
                    TreeMap<Integer, Double> lengthBounds = new TreeMap<>();
                    lengthBounds.put(index, bound);
                    allBounds.put(timeLength, lengthBounds);
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return allBounds;
    }


    public HashMap processesBounds2(String filename) {
        HashMap<Pair, Double> timeBounds = new HashMap<>();
        try {
            Scanner boundScanner = new Scanner(new File(filename));
            while (boundScanner.hasNextLine()) {
                String values[] = boundScanner.nextLine().split(" ");
                int timeLength = Integer.parseInt(values[0]);
                int index = Integer.parseInt(values[1]);
                double bound = Double.parseDouble(values[2]);

                timeBounds.put(new Pair(index, index + timeLength), bound);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return timeBounds;
    }




    public double pruneWithAllBounds(double cond, HashMap<Pair, Double> timeBounds) {
        double prunedBounds = 0;
        double total = 0;
        for (int range = 0; range < 100; range++) {
            for (int n = 0; n < (100 - range); n++) {
                int end = n + range;
                total++;
                if (timeBounds.get(new Pair(n, end)) > cond) prunedBounds++;
            }
        }
        System.out.println("Pruned: " + prunedBounds);
        System.out.println("Total " + total);

        return (prunedBounds / total) * 100;

    }


    public List<Integer> findIntervals(int start, int end, int maxSize) {
        ArrayList<Integer> intervals = new ArrayList<>();
        return findIntervals(intervals, start, end, maxSize);
    }

    public List<Integer> findIntervals(ArrayList<Integer> intervals, int start, int end, int maxSize) {
        int twoPower = maxSize;
        while (twoPower > 0) {
            if (start == end) {
                intervals.add(1);
                return intervals;
            } else if ((start == 0 || start % twoPower == 0) && start + (twoPower - 1) <= end) {
                if (start + (twoPower - 1) == end) {
                    intervals.add(twoPower);
                    return intervals;
                } else {
                    intervals.add(twoPower);
                    return findIntervals(intervals, start + twoPower, end, maxSize);
                }
            }
            twoPower = twoPower == 1 ? -1 : twoPower / 2;
        }

        return intervals;
    }


    public double pruneWithCompositeBounds(Prune prune, double cond, HashMap<Pair, Double> timeBounds) {
        int maxSize = Integer.highestOneBit(prune.totalTime);
        TreeMap<Integer, HashMap<Integer, Double>> nodeGraph = prune.getNodeGraph();

        double prunedBounds = 0;
        double total =  0;

        for (int range = 0; range < 100; range++) {
            for (int n = 0; n < (100 - range); n++) {
                int end = n + range;
                List<Integer> intervals = findIntervals(n, end, maxSize);
                double compositeBound = 0.0;
                int temp = 0;
                for (int interval : intervals) {
                    int subStart = n + temp;
                    int subEnd = subStart + (interval - 1);
                    compositeBound += prune.getIntervalWeight(nodeGraph, subStart, subEnd, n, end) * timeBounds.get(new Pair(subStart, subEnd));
                    temp += interval;
                }

                if (compositeBound > cond) prunedBounds++;
                total++;
            }
        }
        return (prunedBounds/total) * 100;


    }

    public static void main(String[] qwerty) {

        final int RANGE = 9;

        PruneRunner runner = new PruneRunner();

        ArrayList<String> filenames = new ArrayList<>();

//        filenames = runner.addMiscFiles(filenames);
//        filenames = runner.addIGapFiles(filenames, RANGE);
        filenames.add("iGap_series_002_run03.txt");
//        filenames.add("multi_burst0.txt");
//        filenames.add("D4D_January_hourly.txt");
//        filenames.add("simple_example.txt");

        String basepath = "output";
        double average_thresh_conductance = 0d;
        double average_burst_conductance = 0d;
        for (String filename : filenames) {
            Prune prune = runner.loadFile(filename, basepath);
            double burstCond = runner.getMetaInfo(basepath, filename, prune).rawConductance;
            double thresholdedCond = runner.doThresholds(basepath, filename, prune).rawConductance;

            HashMap<Pair, Double> timeBounds = runner.processesBounds2("iGap002run03_bounds.txt");
            System.out.println("Pruning with ALL bounds with BURST: "+runner.pruneWithAllBounds(burstCond, timeBounds));
            System.out.println("Pruning with ALL bounds with THRESH: "+runner.pruneWithAllBounds(thresholdedCond, timeBounds));

            Instant tStart = Instant.now();
            System.out.println("Pruning with COMPOSITE bounds: "+
                runner.pruneWithCompositeBounds(prune, thresholdedCond, timeBounds));

           System.out.println(runner.findIntervals(13, 89, 64));

             Instant tEnd = Instant.now();
             System.out.println("Duration: " + Duration.between(tStart, tEnd).getSeconds());



        }

    }
}
