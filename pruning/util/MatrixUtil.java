package pruning.util;

import pruning.Edge;

import java.util.HashMap;
import java.util.TreeMap;

/**
 * Created by gaurav on 4/24/16.
 */
public class MatrixUtil {

    private TreeMap<Integer, HashMap<Edge, Double>> timeGraph;
    private int totalNodes;

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

    public double getDiagonalWeights(int start, int end) {
        double sum = 0;
        HashMap<Integer, Double> diagonalWeights = new HashMap<>();
        for (int i = start; i <= end; i++) {
            for (Edge e : timeGraph.get(i).keySet()) {
                double weight = timeGraph.get(i).get(e);
                diagonalWeights.compute(e.getNode1(), (k, v) -> v == null ? weight : v + weight);
                diagonalWeights.compute(e.getNode2(), (k, v) -> v == null ? weight : v + weight);
            }
        }
        return diagonalWeights.values().stream().mapToDouble(Double::doubleValue).sum();
    }
}
