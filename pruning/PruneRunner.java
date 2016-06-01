package pruning;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.time.Duration;
import java.time.Instant;

import java.util.*;

public class PruneRunner {

    public Prune loadFile(String filename, String basepath, boolean isOld) {
        return new Prune(filename, basepath, isOld);
    }

    public TreeMap<Integer, List<Set<Edge>>> getConnectedData(Prune prune, double k) {
        return prune.getConnectedComponents(k);
    }

    public Result getResult(Prune prune, TreeMap<Integer, List<Set<Edge>>> connectedComps ) {
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

    public double thresholdBestResult(String filename, String basepath, boolean isOld) {
        Prune prune = loadFile(filename, basepath, isOld);
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

//        iGapFiles.add(fileStart + "_020_run0");
//        iGapFiles.add(fileStart + "_015_run0");
        iGapFiles.add(fileStart + "_010_run0");
//        iGapFiles.add(fileStart + "_005_run0");
//        iGapFiles.add(fileStart + "_002_run0");

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

    public Result doThresholds(Prune prune, double k) {
        System.out.println("For k: " + k);
        TreeMap<Integer, List<Set<Edge>>> connectedComps = getConnectedData(prune, k);
        Result result = getResult(prune, connectedComps);
         prune.printTop5();
        return result;
    }

    public List<Result> doThresholds(Prune prune, double[] kList) {
        List<Result> results = new ArrayList<>();
        for (double k : kList) {
            System.out.println("For k: " + k);
            TreeMap<Integer, List<Set<Edge>>> connectedComps = getConnectedData(prune, k);
            results.add(getResult(prune, connectedComps));
//            prune.printTop5();
        }
        return results;
    }

    public double[][] getAllBounds(Prune prune) {

        int matrixSize = getAllBoundsSize(prune.totalTime);
        double result[][] = new double[matrixSize][3];

        int n;
        int total = 0;
        for (int range = 0; range < prune.totalTime; range++) {
            for (n = prune.startTime; n < (prune.endTime - range); n++) {
                int end = n + range;
                double bounds = prune.getBounds(n, end);
                result[total][0] = n;
                result[total][1] = end;
                result[total][2] = bounds;

//                System.out.println("EXACT "+n + " " + end + " " + bounds);
                total++;

//                System.out.printf(",%f",bounds);
            }
        }
        return result;
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

    public int getAllBoundsSize(int size) {
        int result = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < (size - i); j++) {
                result++;
            }
        }
        return result;
    }


    public double[][] getCompositeBounds(Prune prune) {

        int matrixSize = getAllBoundsSize(prune.totalTime);
        double result[][] = new double[matrixSize][3];

        int totalTime = prune.totalTime;

        // Get the largest power of two that is less or equal to the total time of the whole graph
        int maxSize = Integer.highestOneBit(totalTime);

        TreeMap<Integer, HashMap<Integer, Double>> nodeGraph = prune.getNodeWeights();

        HashMap<Pair, Double> preBounds = new HashMap<>();

        int index = 0;

        for (int range = 0; range < totalTime; range++) {
            for (int n = prune.startTime; n < (prune.endTime - range); n++) {
                int end = n + range;
                List<Integer> intervals = findIntervals(n, end, maxSize);
                double compositeBound = 0.0;
                int temp = 0;
                for (int interval : intervals) {
                    int subStart = n + temp;
                    int subEnd = subStart + (interval - 1);
                    compositeBound += prune.getIntervalWeight(nodeGraph, subStart, subEnd, n, end) *
                        getPrecomputedBounds(prune, subStart, subEnd, preBounds);
                    temp += interval;
                }
                result[index][0] = n;
                result[index][1] = end;
                result[index][2] = compositeBound;

//                System.out.println("COMPOSITE " + n + " " + end + " " + compositeBound);
                System.out.printf("%f\n",compositeBound);

                index++;
            }
        }

        System.out.println("total: "+index);

        return result;
    }

    public double getPrecomputedBounds(Prune prune, int start, int end, HashMap<Pair, Double> preBounds) {
        Pair timePeriod = new Pair(start, end);
        if (preBounds.containsKey(timePeriod)) return preBounds.get(timePeriod);
        else {
            double bounds = prune.getBounds(start, end);
            preBounds.put(timePeriod, bounds);
            return bounds;
        }
    }

    public double getPercentagePruned(double[][] bounds, double cond) {

//        System.out.println("Matrix length: " + bounds.length);
//        System.out.println("Matrix[0] length: " + bounds[0].length);

        double prunedBounds = 0;
        for (int i = 0; i < bounds.length; i++) {
            if (bounds[i][2] > cond) prunedBounds++;
        }
        return (prunedBounds / bounds.length) * 100;
    }

    public static void main(String[] qwerty) {

        PruneRunner runner = new PruneRunner();

        ArrayList<String> filenames = new ArrayList<>();

//        filenames = runner.addMiscFiles(filenames);
//        filenames = runner.addIGapFiles(filenames, RANGE);
//        filenames.add("iGap_series_010_run03.txt");
//        filenames.add("D4D_January_hourly.txt");

//        filenames.add("trade_data.txt");
        filenames.add("mailing_list_data.txt");
//        filenames.add("multiburst_five_varied.txt");


        String basepath = "output";

        for (String filename : filenames) {
            Prune prune = runner.loadFile(filename, basepath, false);

//            double burstCond = runner.getMetaInfo(basepath, filename, prune).rawConductance;
            double thresholdedCond = runner.doThresholds(prune, 6).rawConductance;

//            double[][] allBounds = runner.getAllBounds(prune);

            System.out.println("whole bound: "+prune.getBounds(prune.startTime, prune.endTime));

            Instant tStart = Instant.now();
            double[][] compositeBounds = runner.getCompositeBounds(prune);
            Instant tEnd = Instant.now();
            System.out.println("Duration: " + Duration.between(tStart, tEnd).getSeconds());

//            System.out.println("Percentage pruned with all bounds: "+runner.getPercentagePruned(allBounds, thresholdedCond));
//            System.out.println("Percentage pruned with composite bounds: "+runner.getPercentagePruned(compositeBounds, thresholdedCond));

        }
    }

    /************************************************************************************************************
     public void compareLibs(Prune prune) {

     Instant tStart = Instant.now();
     System.out.println("JBLAS result "+prune.getBounds(15,25));
     System.out.println("JBLAS result "+prune.getBounds(15,35));
     System.out.println("JBLAS result "+prune.getBounds(15,45));
     System.out.println("JBLAS result "+prune.getBounds(15,55));
     System.out.println("JBLAS result "+prune.getBounds(15,65));
     Instant tEnd = Instant.now();
     System.out.println("JBLAS time take "+Duration.between(tStart, tEnd).getSeconds());

     tStart = Instant.now();
     System.out.println("MTJ result "+prune.testMTJ(15,25));
     System.out.println("MTJ result "+prune.testMTJ(15,35));
     System.out.println("MTJ result "+prune.testMTJ(15,45));
     System.out.println("MTJ result "+prune.testMTJ(15,55));
     System.out.println("MTJ result "+prune.testMTJ(15,65));
     tEnd = Instant.now();
     System.out.println("MTJ time take "+Duration.between(tStart, tEnd).getSeconds());
     }
     *************************************************************************************************************/
}
