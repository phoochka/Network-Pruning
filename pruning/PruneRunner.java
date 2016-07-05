package pruning;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.time.Duration;
import java.time.Instant;

import java.util.*;
import java.util.function.DoubleFunction;
import java.util.function.Function;
import java.util.function.IntToDoubleFunction;

public class PruneRunner {

    public Prune loadFile(String filename, Path basepath, boolean isOld) {
        return new Prune(filename, basepath, isOld);
    }

    public TreeMap<Integer, List<Set<Edge>>> getConnectedData(Prune prune, double k) {
        return prune.getConnectedComponents(k);
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

    public double thresholdBestResult(String filename, Path basepath, boolean isOld, IntToDoubleFunction f) {
        Prune prune = loadFile(filename, basepath, isOld);
        final double STARTING_K = 1.0d;
        final double ACCEPTABLE_CONDUCTANCE = 10;
        final double K_INCREMENT = 1;

        double baseResult = prune.traverseGraph(prune.getConnectedComponents(STARTING_K), f).normalizedConductance;
        double incResult = Double.MAX_VALUE;
        double k = STARTING_K;
        while (((baseResult * 0.10) + baseResult) <= incResult) {
            incResult = prune.traverseGraph(prune.getConnectedComponents(k + K_INCREMENT), f).normalizedConductance;
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
                for (int j = 0; j < matrix[i].length; j++) {
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

    public Result doThresholds(Prune prune, double k, IntToDoubleFunction function) {
        System.out.println("For k: "+k);
        TreeMap<Integer, List<Set<Edge>>> connectedComps = getConnectedData(prune, k);
        Result result = prune.traverseGraph(connectedComps, function);
//         prune.printTop(25);
        return result;
    }

    public List<Result> doThresholds(Prune prune, double[] kList, IntToDoubleFunction function) {
        List<Result> results = new ArrayList<>();
        for (double k : kList) {
            System.out.println("For k: " + k);
            TreeMap<Integer, List<Set<Edge>>> connectedComps = getConnectedData(prune, k);
            results.add(prune.traverseGraph(connectedComps,function));
//            prune.printTop5();
        }
        return results;
    }

    public double[][] readBoundsFile(String filename) {
        ArrayList<double[]> list = new ArrayList<>();

        try {
            Scanner boundScanner = new Scanner(new File("bounds/"+filename));
            while(boundScanner.hasNextLine()) {
                String strLine[] = boundScanner.nextLine().split(" ");
                double[] dLine = new double[3];
                for(int i=0;i<3;i++) dLine[i] = Double.parseDouble(strLine[i]);
                list.add(dLine);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        double[][] result = new double[list.size()][3];
        int count = 0;
        for(double[] d : list) {
            result[count][0] = d[0];
            result[count][1] = d[1];
            result[count][2] = d[2];
            count++;
        }

        return result;
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

                total++;

//                System.out.printf("\n%f",bounds);
            }
            System.out.println("All bounds computed range "+range);
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
//                System.out.printf("%f\n",compositeBound);

                index++;
            }
            System.out.println("Composite bounds computed range "+range);
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

    public double getPercentagePruned(double[][] bounds, double cond, IntToDoubleFunction normFunciton) {

//        System.out.println("Matrix length: " + bounds.length);
//        System.out.println("Matrix[0] length: " + bounds[0].length);

        double prunedBounds = 0;
        for (int i = 0; i < bounds.length; i++) {
            int timePeriod = (int) (bounds[i][1] - bounds[i][0] + 1);
            double normalizedBounds = bounds[i][2] * normFunciton.applyAsDouble(timePeriod);
            if (normalizedBounds > cond) prunedBounds++;
        }
        return (prunedBounds / bounds.length) * 100;
    }

    public double checkRootSolution(Prune prune, double alpha, double k, String filename) {

        double root = 1.0 / alpha;
        IntToDoubleFunction function = (i) -> 1/Math.pow(i, root);
        Instant tStart = Instant.now();
        double solution = doThresholds(prune, k, function).normalizedConductance;
        Instant tEnd = Instant.now();
        System.out.println("Duration: " + Duration.between(tStart, tEnd).getSeconds());
        double[][] allBounds = readBoundsFile("results/all_bounds_"+filename);

        System.out.println("Percentage pruned for: "+filename+" with alpha "+alpha+" and k "+k+" = "+getPercentagePruned(allBounds, solution, function));

        return solution;
    }
    public double checkSolution(Prune prune, double alpha, double k, String filename) {
        IntToDoubleFunction function = (i) -> Math.exp(-alpha * i);
        Instant tStart = Instant.now();
        double solution = doThresholds(prune, k, function).normalizedConductance;
        Instant tEnd = Instant.now();
        System.out.println("Duration: " + Duration.between(tStart, tEnd).getSeconds());
        double[][] allBounds = readBoundsFile("all_bounds_"+filename);

        System.out.println("Percentage pruned for: "+filename+" with alpha "+alpha+" and k "+k+" = "+getPercentagePruned(allBounds, solution, function));

        return solution;
    }

    public static void main(String[] qwerty) {


        PruneRunner runner = new PruneRunner();

        runner.findIntervals(11, 31, 64);

        System.exit(0);
        Path basepath = FileSystems.getDefault().getPath("output");

        ArrayList<String> filenames = new ArrayList<>();

//        filenames = runner.addMiscFiles(filenames);
//        filenames = runner.addIGapFiles(filenames, RANGE);
//        filenames.add("iGap_series_020_run03.txt");
//        filenames.add("iGap_series_010_run05.txt");
        filenames.add("iGap_series_002_run03.txt");

//        filenames.add("trade_data.txt");
//        filenames.add("mailing_list_data.txt");
//        filenames.add("trade_data_log.txt");

//        filenames.add("network_traffic_GCC.txt");
//        filenames.add("GDELT_GCC_2010_weekly.txt");
//        filenames.add("GDELT_GCC_2011_weekly.txt");


//        filenames.add("intWeightDifficult2_0500_run004.txt");

        for (String filename : filenames) {
            Prune prune = runner.loadFile(filename, basepath, true);

//            double burstCond = runner.getMetaInfo(basepath, filename, prune).rawConductance;

            IntToDoubleFunction cub = i -> 1.0/Math.cbrt(i);

            runner.printResults(prune, filename);


//            runner.checkSolution(prune, 0.01, 10, filename);
//            runner.checkSolution(prune, 0.02, 10, filename);
//            runner.checkSolution(prune, 0.04, 10, filename);
//            runner.checkSolution(prune, 0.08, 10, filename);
//            runner.checkSolution(prune, 0.10, 10, filename);
//            runner.checkSolution(prune, 0.12, 10, filename);
//            runner.checkSolution(prune, 0.15, 10, filename);
//

//            Instant tStart = Instant.now();
//            double[][] allBounds = runner.getAllBounds(prune);
//            Instant tEnd = Instant.now();
//            System.out.println("Duration to find all bounds for "+filename+": " + Duration.between(tStart, tEnd).getSeconds());
//
//
//            tStart = Instant.now();
//            double[][] compositeBounds = runner.getCompositeBounds(prune);
//            tEnd = Instant.now();
//            System.out.println("Duration to find composite bounds for "+filename+": " + Duration.between(tStart, tEnd).getSeconds());
//
//            runner.printMatrix("bounds/all_bounds_"+filename,allBounds);
//            runner.printMatrix("bounds/com_bounds_"+filename,compositeBounds);
//
//            double cubeSolution = runner.doThresholds(prune, 10, cub).normalizedConductance;
//            System.out.println("Percentage pruned with all bounds with cube solution: "+runner.getPercentagePruned(runner.readBoundsFile("all_bounds_"+filename), cubeSolution, cub));

//            IntToDoubleFunction function = (i) -> Math.exp(-0.006 * i);
//            double hashingSolution1 = 0.0004561284;
//            System.out.println("Percentage pruned with all bounds with hashing" +
//                    " solution: "+runner.getPercentagePruned(runner.readBoundsFile("all_bounds_"+filename), hashingSolution1, cub));

//            double hashingSolution2 = 0.1388448605;
//            System.out.println("Percentage pruned with all bounds with hashing solution: "+runner.getPercentagePruned(runner.readBoundsFile("all_bounds_"+filename), hashingSolution2, cub));
//            System.out.println("Percentage pruned with composite bounds: "+runner.getPercentagePruned(compositeBounds, t1));

        }
    }

    private void printResults(Prune prune, String filename) {
       try {
        PrintWriter writer = new PrintWriter("bounds/root_solution_"+filename);
            for (double i = 1.4; i < 5.0; i+= 1.0){
                writer.write(i+"\t"+checkRootSolution(prune, i, 10, filename)+"\n");
            }
           writer.close();
       } catch (FileNotFoundException e) {
           e.printStackTrace();
       }
    }
}
