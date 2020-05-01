import java.util.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.time.Instant;
import java.time.Duration;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JButton;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.Box;
import javax.swing.JOptionPane;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.geom.Line2D;

public class main {
    private static Population currentPopulation;
    private static Population nextPopulation;
    private static GUI gui;
    private static int numOfGenerations;
    private static int currentGeneration;

    public static void main(String[] args) {
        ArrayList<String> edgeStrings = getNodesFromFile("input.txt");
        NodeList nodes = generateNodeList(edgeStrings);
        System.out.println("Graph Adjacency Matrix");
        System.out.println();
        nodes.printAdjacencyMatrix();
        System.out.println();
        System.out.println();

        numOfGenerations = getInput("Number of generations", 10000);
        int populationSize = getInput("Population Size", 1000);
        int crossoverRate = getInput("Crossover Rate", 100);
        int mutationRate = getInput(
                "Mutation Rate (less than " + (100 - crossoverRate) + ")",
                100 - crossoverRate);


        //System.out.println("Generating random orderings");
        currentPopulation = new Population(nodes.size(), populationSize, crossoverRate,
                mutationRate, nodes);
        currentPopulation.generateRandomOrderings();
        gui = new GUI(currentPopulation, nodes, 5, new ButtonListener());
    }

    static class ButtonListener implements ActionListener {
        public ButtonListener(){};

        @Override
        public void actionPerformed(ActionEvent e) {
            if(e.getActionCommand().equals("Next Generation")) {
                if (currentPopulation.getBestFromGeneration().isSolved()) {
                    System.out.println("Solved at generation: " + currentGeneration);
                } else {
                    nextGeneration();
                    Ordering o = currentPopulation.getBestFromGeneration();
                    gui.update(o, o.getFitnessCost(), ++currentGeneration);

                }
            } else if(e.getActionCommand().equals("Last Generation")) {
                for(;currentGeneration < numOfGenerations; ++currentGeneration) {
                    nextGeneration();
                    Ordering o = currentPopulation.getBestFromGeneration();
                    if(o.isSolved()) {
                        System.out.println("Solved at generation: " + currentGeneration);
                        break;
                    }
                }
                Ordering o = currentPopulation.getBestFromGeneration();
                gui.update(o, o.getFitnessCost(), currentGeneration);
            } else if(e.getActionCommand().equals("Reset")) {
                reset();
            } else if(e.getActionCommand().equals("Swap Fitness Function")) {
                reset();
                if(currentPopulation.getFitnessFunction() == Population.FitnessFunction.TETTAMANZI) {
                    currentPopulation.setFitnessFunction(Population.FitnessFunction.DEFAULT);
                    gui.setFitnessFunctionText(Population.FitnessFunction.DEFAULT);
                } else {
                    currentPopulation.setFitnessFunction(Population.FitnessFunction.TETTAMANZI);
                    gui.setFitnessFunctionText(Population.FitnessFunction.TETTAMANZI);
                }
            }
        }
    }

    private static void reset() {
        currentPopulation.generateRandomOrderings();
        nextPopulation = null;
        currentGeneration = 0;
        Ordering startOrdering = currentPopulation.getOrdering(0);
        gui.update(startOrdering, startOrdering.getFitnessCost(), currentGeneration);
    }

    private static void nextGeneration() {
        if(nextPopulation == null) {
            nextPopulation = currentPopulation.generateNextPopulation();
        } else {
            currentPopulation = nextPopulation;
            nextPopulation = currentPopulation.generateNextPopulation();
        }
    }

    private static ArrayList<String> getNodesFromFile(String path) {
        ArrayList<String> r;
        try {
            r = new ArrayList(Files.readAllLines(Paths.get(path)));
        } catch (IOException e) {
            return null;
        }
        return r;
    }

    /*
     * generateNodeList takes the string from the input file in the form 12 18
     * where the first value is the head of the edge and the second is the tail
     */
    private static NodeList generateNodeList(ArrayList<String> edges) {
        NodeList nodes = new NodeList();
        for (String edgeNodes : edges) {
            String[] headAndTail = edgeNodes.split(" ");
            Integer headNum = Integer.parseInt(headAndTail[0]);
            Integer tailNum = Integer.parseInt(headAndTail[1]);

            Node head = new Node(headNum);
            Node tail = new Node(tailNum);
            // Check if head of edge exists
            if(nodes.contains(head)) {
                head = nodes.getNode(head);
                head.addTail(tailNum);
            } else {
                nodes.add(head);
                head.addTail(tailNum);
            }

            if(nodes.contains(tail)) {
                tail = nodes.getNode(tail);
                tail.addTail(headNum);
            } else {
                nodes.add(tail);
                tail.addTail(headNum);
            }
        }
        return nodes;
    }

    private static int getInput(String inputTitle, int max) {
        String inputString = JOptionPane.showInputDialog(inputTitle);
        if (inputString == null) System.exit(0);
        while (!isInputValid(inputString, max)) {
            inputString = JOptionPane.showInputDialog(inputTitle);
            if (inputString == null) System.exit(0);
        }
        return Integer.parseInt(inputString);
    }

    private static boolean isInputValid(String s, int max) {
        try {
            int i = Integer.parseInt(s);
            if (0 < i && i < max) {
                return true;
            }
            throw new NumberFormatException();
        } catch(NumberFormatException e) {
            JOptionPane.showMessageDialog(null,
                    "Value must be an integer less than " + max);
        }
        return false;
    }
}


/**
 * Orderings class contains an ordering of nodes
 * Each ordering has an order which describes the order of nodes in the ordering,
 * a length which describes the length of the ordering
 * and a fitness cost which is the fitness cost of the ordering
 */
class Ordering {
    private int[] order;
    public int length;
    private double fitnessCost;

    public Ordering() {
        order = null;
    }

    public Ordering(ArrayList<Integer> ordering) {
        length = ordering.size();
        order = new int[ordering.size()];
        for(int i = 0; i < order.length; i++) {
            order[i] = ordering.get(i);
        }
    }

    public Ordering(int[] ordering) {
        this.order = ordering;
        length = ordering.length;
    }

    public int[] getOrderingArray() {
        return order;
    }

    public int get(int index) {
        if (index >= 0 && index < length) {
            return order[index];
        }
        // TODO throw exception
        return -1;
    }


    public void setFitnessCost(double fitnessCost) {
        this.fitnessCost = fitnessCost;
    }

    public double getFitnessCost() {
        return fitnessCost;
    }

    // Point is represented as a double array where x = double[0], y=double[1]
    // Generates map from node value to node coordinate
    public Map<Integer, double[]> generateCoordinateMap(double chunkValue) {
        Map<Integer, double[]> coordinates = new HashMap<Integer, double[]>();
        for(int i = 0; i < length; i++) {
            double x = Math.cos(i * chunkValue);
            double y = Math.sin(i * chunkValue);
            coordinates.put(order[i], new double[]{x, y});
        }
        return coordinates;
    }

    /**
     * removeDuplicates removes duplicates from an ordering and replaces the duplicates
     * with the values that are missing in the ordering
     */
    public void removeDuplicates() {
        int[] rollCall = new int[order.length];
        for (int i : order) {
            rollCall[i] += 1;
        }

        for (int nodeIndex = 0; nodeIndex < order.length; nodeIndex++) {
            if (rollCall[order[nodeIndex]] > 1) {
                int replacement;
                for (replacement = 0; replacement < rollCall.length; replacement++) {
                    if (rollCall[replacement] == 0) {
                        rollCall[replacement] += 1;
                        break;
                    }
                }

                // Lower the value in roll call to 1 before replacement otherwise replacement
                // will be decremented again and will still be counted as missing
                rollCall[order[nodeIndex]] -= 1;
                order[nodeIndex] = replacement;
            }
        }
    }

    public boolean isSolved() {
        Integer[] solvedOrder = {5, 0, 6, 15, 4, 13, 1, 11, 10, 14, 12, 9, 17, 7, 3, 2, 16, 8};
        Integer[] solvedOrder2 = {5, 6, 0, 15, 4, 13, 1, 11, 10, 14, 12, 9, 17, 7, 3, 2, 16, 8};
        Integer[] orderAsIntegers = Arrays.stream(order).boxed().toArray(Integer[]::new);

        int startIndex = -1;
        boolean solved = true;
        boolean started = false;
        for(int i = 0; i < orderAsIntegers.length; i++) {
            if(started) break;
            if(orderAsIntegers[i] == solvedOrder[0]) {
                startIndex = i; // If not found going forward we can skip initial search when looking back
                started = true;
                for(int j = 1; j < solvedOrder.length; j++) {
                    if(orderAsIntegers[(j + i) % solvedOrder.length] != solvedOrder[j]
                            && orderAsIntegers[(j + i) % solvedOrder.length] != solvedOrder2[j]) {
                        solved = false;
                        break;
                    }
                }
            }
        }

        if(!solved) {
            solved = true;
            for(int j = 1; j < solvedOrder.length; j++) {
                if(startIndex - j < 0) startIndex = startIndex + orderAsIntegers.length;
                if(orderAsIntegers[startIndex-j] != solvedOrder[j] && orderAsIntegers[startIndex-j] != solvedOrder2[j]) {
                    solved = false;
                    break;
                }
            }
        }
        return solved;
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        Ordering clone = new Ordering(order.clone());
        clone.setFitnessCost(fitnessCost);
        return clone;
    }

    @Override
    public String toString() {
        String output = "";
        output += "{ ";
        for(int i : order) output += (i + " ");
        return output += "}";
    }
}


/**
 * Population class contains a generation population
 * This class generates the next population by using
 * selection, reproduction, mutation and crossover
 */
class Population {
    private int graphSize;
    private int populationSize;
    private int crossoverRate;
    private int mutationRate;
    private double chunk;
    private NodeList nodelist;
    private Ordering[] orderings;

    private static final double INTERSECTION_PARAM = 0.9;
    private static final double LENGTH_DEVIATION_PARAM = 0.85;
    private static final double ANGLE_DEVIATION_PARAM = 0.20;

    public static enum FitnessFunction {
        DEFAULT,
        TETTAMANZI
    }

    private FitnessFunction fitnessFunction = FitnessFunction.DEFAULT;

    public Population(int graphSize, int populationSize,
            int crossoverRate, int mutationRate, NodeList nodelist) {
        this.graphSize = graphSize;
        this.populationSize = populationSize;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.nodelist = nodelist;
        orderings = new Ordering[populationSize];
        chunk = (2 * Math.PI) / graphSize;
    }

    public Population(int graphSize, int populationSize,
            int crossoverRate, int mutationRate, NodeList nodelist, Population.FitnessFunction fitnessFunc) {
        this.graphSize = graphSize;
        this.populationSize = populationSize;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.nodelist = nodelist;
        orderings = new Ordering[populationSize];
        chunk = (2 * Math.PI) / graphSize;
        fitnessFunction = fitnessFunc;
    }

    public void setOrderings(Ordering[] orderings) {
        this.orderings = orderings;
    }

    public void setFitnessFunction(FitnessFunction ff) {
        fitnessFunction = ff;
    }

    public FitnessFunction getFitnessFunction() {
        return fitnessFunction;
    }

    public Ordering getOrdering(int index) {
        if(index < orderings.length) {
            return orderings[index];
        }
        return null;
    }

    public double getChunk() {
        return chunk;
    }

    /**
     * generateRandomOrderings is the initial method used during generation 0 to
     * create the initial Orderings
     */
    public void generateRandomOrderings() {
        ArrayList<Integer> baseNumbers = new ArrayList<>();
        for(int i = 0; i < graphSize; i++) {
            baseNumbers.add(i);
        }

        for(int i = 0; i < orderings.length; i++) {
            Collections.shuffle(baseNumbers);
            orderings[i] = new Ordering(baseNumbers);
            System.out.println((i + 1)  + " => " + orderings[i]);
        }
    }

    /**
     * generateNextPopulation method responsible for creating the next population from
     * the current one. Calculates fitness cost for all orderings, performs selection
     * algorithm on current orderings and finally iterates through the current orderings
     * until there arent any left to use and chooses to either mutate, crossover or
     * reproduce from the selected orderings and adds the output of the algorithms
     * to a new generation of orderings.
     */
    public Population generateNextPopulation() {
        calculateFitnessForAllOrderings();
        if(fitnessFunction == FitnessFunction.TETTAMANZI) {
            performSelectionT();
        } else {
            performSelection();
        }

        // Generate the next generation
        Ordering[] nextGenerationOrderings = new Ordering[populationSize];
        // Records which orderings from the previous population have been used
        ArrayList<Integer> orderPositionsRemaining = new ArrayList<>();
        for(int i = 0; i < populationSize; i++) { orderPositionsRemaining.add(i);}
        Random rand = new Random();

        int fillOrderingsIndex = 0;
        while(orderPositionsRemaining.size() > 0) {
            int probability = rand.nextInt(101);
            if(crossoverRate >= probability) {
                // Unable to perform crossover unless there are at least two
                // orderings left to choose from
                if(orderPositionsRemaining.size() > 1) {
                    int orderRemainingIndex1 = rand.nextInt(
                            orderPositionsRemaining.size());
                    int orderRemainingIndex2 = rand.nextInt(
                            orderPositionsRemaining.size());
                    while(orderRemainingIndex1 == orderRemainingIndex2) {
                        orderRemainingIndex2 = rand.nextInt(orderPositionsRemaining.size());
                    }

                    int orderIndex1 = orderPositionsRemaining.get(orderRemainingIndex1);
                    int orderIndex2 = orderPositionsRemaining.get(orderRemainingIndex2);
                    ArrayList<Ordering> crossoverOrderings
                        = crossover(orderIndex1, orderIndex2);

                    for(Ordering ordering : crossoverOrderings) {
                        nextGenerationOrderings[fillOrderingsIndex] = ordering;
                        fillOrderingsIndex++;
                    }

                    // Removing two indexes can cause an out of bounds error
                    // Solution is to remove by largest index first
                    if (orderRemainingIndex1 > orderRemainingIndex2) {
                        orderPositionsRemaining.remove(orderRemainingIndex1);
                        orderPositionsRemaining.remove(orderRemainingIndex2);
                    } else {
                        orderPositionsRemaining.remove(orderRemainingIndex2);
                        orderPositionsRemaining.remove(orderRemainingIndex1);
                    }

                }
            } else if (probability < (crossoverRate + mutationRate)) {
                int orderRemainingIndex = rand.nextInt(orderPositionsRemaining.size());
                int orderIndex = orderPositionsRemaining.get(orderRemainingIndex);
                nextGenerationOrderings[fillOrderingsIndex] = mutate(orderIndex);
                orderPositionsRemaining.remove(orderRemainingIndex);
                fillOrderingsIndex++;
            } else if((crossoverRate + mutationRate) <= probability) {
                int orderRemainingIndex = rand.nextInt(orderPositionsRemaining.size());
                int orderIndex = orderPositionsRemaining.get(orderRemainingIndex);
                nextGenerationOrderings[fillOrderingsIndex] = reproduce(orderIndex);
                orderPositionsRemaining.remove(orderRemainingIndex);
                fillOrderingsIndex++;
            }
        }

        Population nextPopulation = new Population(graphSize, populationSize, crossoverRate,
                mutationRate, nodelist, fitnessFunction);
        nextPopulation.setOrderings(nextGenerationOrderings);

        //System.out.println("Next Gen orderings");
        //for(int i = 0; i < nextGenerationOrderings.length; i++) {
            //System.out.println((i + 1)  + " => " + nextGenerationOrderings[i]);
        //}

        return nextPopulation;
    }

    /**
     * calculateFitnessForAllOrderings loops through all orderings in the orderings
     * array and calculates that orderings fitness cost
     */
    private void calculateFitnessForAllOrderings() {
        long totalTime = 0;
        if(fitnessFunction == FitnessFunction.TETTAMANZI) {
            for(Ordering ordering : orderings) {
                totalTime += calculateOrderingFitness2(ordering);
            }
            double avgTimeInMillis = (totalTime / orderings.length) / 100000.0;
            //System.out.println("Average time to compute tettamanzi fitness: " + avgTimeInMillis + "ms");
        } else {
            for(Ordering ordering : orderings) {
                totalTime += calculateOrderingFitness(ordering, chunk);
            }
            double avgTimeInMillis = (totalTime / orderings.length) / 100000.0;
            //System.out.println("Average time to compute fitness: " + avgTimeInMillis + "ms");
        }
    }

    /**
     * performSelection removes poor orderings (high fitness cost) by replacing the
     * worst third of the orderings with a copy of the first third of orderings
     */
    private void performSelection() {
        Arrays.sort(orderings, new Comparator<Ordering>() {
            @Override
            public int compare(Ordering o1, Ordering o2) {
                return Double.compare(o1.getFitnessCost(), o2.getFitnessCost());
            }
        });

        Ordering[] topThird = Arrays.copyOfRange(orderings, 0, populationSize / 3);
        int lastThirdStartIndex = orderings.length - topThird.length;
        for(int i = lastThirdStartIndex, j = 0; i < populationSize; i++, j++) {
            orderings[i] = topThird[j];
        }
    }

    /**
     * performSelection removes poor orderings (high fitness cost) by replacing the
     * worst third of the orderings with a copy of the first third of orderings
     */
    private void performSelectionT() {
        Arrays.sort(orderings, new Comparator<Ordering>() {
            @Override
            public int compare(Ordering o1, Ordering o2) {
                if(o1.getFitnessCost() < o2.getFitnessCost()) {
                    return 1;
                } else if(o2.getFitnessCost() < o1.getFitnessCost()) {
                    return -1;
                } else {
                    return 0;
                }
            }
        });

        Ordering[] topThird = Arrays.copyOfRange(orderings, 0, populationSize / 3);
        int lastThirdStartIndex = orderings.length - topThird.length;
        for(int i = lastThirdStartIndex, j = 0; i < populationSize; i++, j++) {
            orderings[i] = topThird[j];
        }
    }

    /**
     * CalculateOrderingFitness calculates the fitness for a single ordering
     * by summing the length of all edges
     */
    private long calculateOrderingFitness(Ordering ordering, double chunk) {
        Instant start = Instant.now();
        double fitnessCost = 0;
        Map<Integer, double[]> coordinates = ordering.generateCoordinateMap(chunk);

        for(int i = 0; i < ordering.length; i++) {
            Node node = nodelist.getNode(ordering.get(i));
            ArrayList<Integer> nodeConnections = node.getTails();

            for(Integer nodeValue : nodeConnections) {
                // to ensure connection distances are only counted once dont calculate distance of nodes less than current node value
                // as they are assumed to already have been calculated
                if(node.getNum() < nodeValue) {
                    Node connectedNode = nodelist.getNode(nodeValue);
                    double[] nodePoint = coordinates.get(node.getNum());
                    double[] connectionPoint = coordinates.get(connectedNode.getNum());
                    double distanceBetweenNodes = distance(nodePoint[0], nodePoint[1],
                            connectionPoint[0], connectionPoint[1]);
                    fitnessCost += distanceBetweenNodes;
                }
            }
        }

        Instant end = Instant.now();
        long timeDiff = Duration.between(start, end).toNanos();
        ordering.setFitnessCost(fitnessCost);
        return timeDiff;
    }

    /**
     * distance calculates the pythagorean distance between two points,
     * used to calculate the fitness cost of an orderings
     */
    public double distance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1,2) + Math.pow(y2 - y1,2));
    }

    /**
     * crossover applies the crossover algorithm to two orderings in the ordering
     * array and returns two new orderings created from the crossover in an ArrayList
     */
    public ArrayList<Ordering> crossover(int orderIndex1, int orderIndex2) {

        Ordering ordering1 = orderings[orderIndex1];
        Ordering ordering2 = orderings[orderIndex2];
        int firstItem = 1;
        int seconLastItem = ordering1.length - 2;
        Random rand = new Random();
        // cutting point from [1, N - 2]
        int cuttingPoint = rand.nextInt(seconLastItem) + firstItem;

        Ordering[] splitOrdering1 = splitOrdering(ordering1, cuttingPoint);
        Ordering[] splitOrdering2 = splitOrdering(ordering2, cuttingPoint);

        Ordering mixedOrdering1 = mergeOrderings(splitOrdering1[0], splitOrdering2[1]);
        Ordering mixedOrdering2 = mergeOrderings(splitOrdering2[0], splitOrdering1[1]);
        mixedOrdering1.removeDuplicates();
        mixedOrdering2.removeDuplicates();

        ArrayList<Ordering> orderingList = new ArrayList<Ordering>();
        orderingList.add(mixedOrdering1);
        orderingList.add(mixedOrdering2);
        return orderingList;
    }

    /**
     * Mutate performs mutation on an ordering in a population
     */
    public Ordering mutate(int orderingIndex) {
        Ordering ordering = orderings[orderingIndex];
        Random rand = new Random();
        int swapIndex1 = rand.nextInt(ordering.length);
        int swapIndex2 = rand.nextInt(ordering.length);
        Ordering mutatedOrdering = mutateOrdering(ordering, swapIndex1, swapIndex2);
        return mutatedOrdering;
    }

    /**
     * reproduce performs reproduction on an ordering in a population
     */
    public Ordering reproduce(int orderingIndex) {
        Ordering ordering = null;
        try {
            ordering = (Ordering)orderings[orderingIndex].clone();
        } catch(CloneNotSupportedException e) {
            e.printStackTrace();
        }
        return ordering;
    }

    /**
     * mutateOrdering takes an ordering and swaps two elements positions in
     * the orderings order array
     */
    private Ordering mutateOrdering(Ordering ordering, int index1, int index2) {
        int[] orderingArray = ordering.getOrderingArray().clone();
        int temp = orderingArray[index1];
        orderingArray[index1] = orderingArray[index2];
        orderingArray[index2] = temp;
        return new Ordering(orderingArray);
    }

    /**
     * splitOrdering cuts an ordering in two at the cutting point and returns
     * the split ordering and two sub-orderings in an array
     */
    private Ordering[] splitOrdering(Ordering ordering, int cuttingPoint) {
        Ordering[] subOrderings = new Ordering[2];
        int[] orderingArray = ordering.getOrderingArray();

        int[] subOrderingArray1 = Arrays.copyOfRange(orderingArray, 0, cuttingPoint);
        int[] subOrderingArray2 = Arrays.copyOfRange(orderingArray, cuttingPoint,
                orderingArray.length);
        subOrderings[0] = new Ordering(subOrderingArray1);
        subOrderings[1] = new Ordering(subOrderingArray2);
        return subOrderings;
    }

    private Ordering mergeOrderings(Ordering ordering1, Ordering ordering2) {
        int[] orderingArray1 = ordering1.getOrderingArray();
        int[] orderingArray2 = ordering2.getOrderingArray();

        ArrayList<Integer> mergedOrdering = new ArrayList<>();
        for(int i : orderingArray1) {
            mergedOrdering.add(i);
        }

        for(int i : orderingArray2) {
            mergedOrdering.add(i);
        }
        return new Ordering(mergedOrdering);
    }

    public Ordering getBestFromGeneration() {
        Ordering[] orderingsCopy = new Ordering[populationSize];
        for (int i = 0; i < orderings.length; i++) {
            Ordering ordering = null;
            try {
                ordering = (Ordering)orderings[i].clone();
                orderingsCopy[i] = ordering;
            } catch(CloneNotSupportedException e) {
                e.printStackTrace();
            }
        }

        Arrays.sort(orderingsCopy, new Comparator<Ordering>() {
            @Override
            public int compare(Ordering o1, Ordering o2) {
                if(fitnessFunction == FitnessFunction.DEFAULT) {
                    return Double.compare(o1.getFitnessCost(), o2.getFitnessCost());
                } else {
                    if(o1.getFitnessCost() < o2.getFitnessCost()) {
                        return 1;
                    } else if(o2.getFitnessCost() < o1.getFitnessCost()) {
                        return -1;
                    } else {
                        return 0;
                    }
                }
            }
        });
        return orderingsCopy[0];
    }

    /* Below here are all the methods needed to implement
     * the second fitness function based on the Tettamanzi
     */
    public long calculateOrderingFitness2(Ordering ordering){
        Instant start = Instant.now();
        int crossings = countEdgeCrossings(ordering);
        double meanRelativeSquareError = calculateMeanRelativeSquareError(ordering);
        double cumulativeSquareAngleDeviation = calculateCumulativeSquareDeviation(ordering);
        double orderingFitness = (LENGTH_DEVIATION_PARAM*(1.0 / (meanRelativeSquareError + 1.0)))
            + (INTERSECTION_PARAM*(1.0 / (crossings + 1.0)))
            + (ANGLE_DEVIATION_PARAM*(1.0 / (cumulativeSquareAngleDeviation + 1.0)));
        Instant end = Instant.now();
        long timeDiff = Duration.between(start, end).toNanos();

        ordering.setFitnessCost(orderingFitness);
        return timeDiff;
    }

    public int countEdgeCrossings(Ordering ordering){
        int intersections = 0;
        for (int i = 0; i < ordering.length; i++) {
            Node node = nodelist.getNode(i);

            for(int tailVal : node.getTails()) {
                if(tailVal < node.getNum()) continue;
                intersections += countEdgeIntersections(node.getNum(), tailVal, ordering);
            }
        }
        // Dividing by two because each intersection is counted twice, this is the simplest way
        // of removing double counting
        return intersections / 2;
    }

    private int countEdgeIntersections(int lineStart, int lineEnd, Ordering o) {
        int intersections = 0;
        Map<Integer, double[]> coordinates = o.generateCoordinateMap(chunk);
        double[] lineStartPoint = coordinates.get(lineStart);
        double[] lineEndPoint = coordinates.get(lineEnd);
        for (int i = 0; i < o.length; i++) {
            if(i == lineStart || i == lineEnd) continue;

            Node node = nodelist.getNode(i);
            for(int tail : node.getTails()) {
                if(tail < node.getNum()) continue;
                if(tail == lineStart || tail == lineEnd) continue;

                double[] comparisonPointStart = coordinates.get(node.getNum());
                double[] comparisonPointEnd = coordinates.get(tail);
                if (Line2D.linesIntersect(
                            lineStartPoint[0], lineStartPoint[1], lineEndPoint[0],
                            lineEndPoint[1], comparisonPointStart[0], comparisonPointStart[1],
                            comparisonPointEnd[0], comparisonPointEnd[1])) {
                    intersections++;
                            }
            }
        }
        return intersections;
    }
    public double calculateMeanRelativeSquareError(Ordering ordering){
        double error = 0;

        Map<Integer, double[]> coordinates = ordering.generateCoordinateMap(chunk);

        for(int i = 0; i < ordering.length; i++) {
            Node node = nodelist.getNode(ordering.get(i));
            ArrayList<Integer> nodeConnections = node.getTails();

            for(Integer nodeValue : nodeConnections) {
                // to ensure connection distances are only counted once dont calculate distance of nodes less than current node value
                // as they are assumed to already have been calculated
                if(node.getNum() < nodeValue) {
                    Node connectedNode = nodelist.getNode(nodeValue);
                    double[] nodePoint = coordinates.get(node.getNum());
                    double[] connectionPoint = coordinates.get(connectedNode.getNum());
                    double distanceBetweenNodes = distance(nodePoint[0], nodePoint[1],
                            connectionPoint[0], connectionPoint[1]);

                    error += Math.pow(((distanceBetweenNodes - 0.347296) / 0.347296), 2);
                }
            }
        }
        return error;
    }

    public double calculateCumulativeSquareDeviation(Ordering ordering){
        double deviation = 0;

        Map<Integer, double[]> coordinates = ordering.generateCoordinateMap(chunk);

        for(int i = 0; i < ordering.length; i++) {
            Node node = nodelist.getNode(i);
            for(int orderNum = 0; orderNum < ordering.length; orderNum++) {
                int firstNode = ordering.get(orderNum);
                if(node.getTails().indexOf(firstNode) == -1) continue;
                //System.out.println("First node found: " + firstNode);

                for(int nextNode = orderNum + 1; nextNode < ordering.length; nextNode++) {
                    int secondNode = ordering.get(nextNode);
                    if(node.getTails().indexOf(secondNode) == -1) continue;
                    double[] translatedFirstNode = translatePoint(
                            coordinates.get(firstNode), coordinates.get(node.getNum()));
                    double[] translatedSecondNode = translatePoint(
                            coordinates.get(secondNode), coordinates.get(node.getNum()));
                    double angle = angleBetweenTwoEdges(translatedFirstNode, translatedSecondNode);
                    deviation += Math.pow(angle - (2*Math.PI) / node.getTails().size(), 2);
                    break;
                }
            }
        }
        return deviation;
    }

    private double angleBetweenTwoEdges(double[] edge1, double[] edge2) {
        double edge1Magnitude = Math.sqrt(Math.pow(edge1[0], 2) + Math.pow(edge1[1], 2));
        double edge2Magnitude = Math.sqrt(Math.pow(edge2[0], 2) + Math.pow(edge2[1], 2));
        double dotProduct = edge1[0]*edge2[0] + edge1[1]*edge2[1];

        double angle = Math.toDegrees(
                Math.acos(dotProduct / (edge1Magnitude * edge2Magnitude)));
        return angle;
    }

    private double[] translatePoint(double[] point, double[] offset) {
        return new double[]{point[0] - offset[0], point[1] - offset[1]};
    }

}


/**
 * NodeList class is a wrapper around an ArrayList that contains Nodes
 * This class provides extra functionality on top of the Node ArrayList
 */
class NodeList {
    private ArrayList<Node> nodes;

    public NodeList() {
        nodes = new ArrayList<>();
    }

    public int size() {
        return nodes.size();
    }

    public void add(Node n) {
        nodes.add(n);
    }

    public void addNewNode(Integer num) {
        Node n = new Node(num);
        nodes.add(n);
    }

    public boolean contains(Node n) {
        return nodes.contains(n);
    }

    public ArrayList<Node> getNodes() {
        return nodes;
    }

    public Node getNode(Integer num) {
        for (Node node : nodes) {
            if(node.getNum() == num) {
                return node;
            }
        }
        return null;
    }

    public Node getNode(Node n) {
        for (Node node : nodes) {
            if(node.getNum() == n.getNum()) {
                return node;
            }
        }
        return null;
    }

    public void printAdjacencyMatrix() {
        HashMap<Integer, ArrayList<Integer>> nodeTailMap = new HashMap<Integer, ArrayList<Integer>>();
        for(Node n : nodes) {
            nodeTailMap.put(n.getNum(), n.getTails());
        }

        //Print top bar of adjacency matrix
        System.out.printf("%3s", " |");
        //nodeTailMap.forEach((k, v) -> System.out.printf("%3s", k + "|"));
        for (int i = 0; i < nodes.size(); i++) System.out.printf("%3s", i + "|");
        System.out.println();
        System.out.println("--+".repeat(19));


        for(int i = 0; i < nodes.size(); i++) {
            System.out.printf("%3s", i + "|");
            for(int j = 0; j < nodes.size(); j++) {
                int connection = (nodeTailMap.get(i).contains(j)) ? 1 : 0;
                System.out.printf("%3s", connection + "|");
            }
            System.out.println();
            System.out.println("--+".repeat(19));
        }
    }
}


/**
 * Node is a class that represents a node from the input file
 * Each node has a number that identifies it and a list of numbers that
 * corresponds to the node values of the nodes it is connected to.
 */
class Node {
    private int num;
    private ArrayList<Integer> tails;

    public Node(int num) {
        this.num = num;
        this.tails = new ArrayList<>();
    }

    public Node(String num) {
        this.num = Integer.parseInt(num);
        this.tails = new ArrayList<>();
    }

    public int getNum() {
        return num;
    }

    public void addTail(Integer num) {
        tails.add(num);
    }

    public ArrayList<Integer> getTails() {
        return tails;
    }

    @Override
    public String toString() {
        return this.num + " " + this.tails.toString();
    }

    @Override
    public boolean equals(Object o) {
        if(o == this) return true;
        if(!(o instanceof Node)) return false;

        Node n = (Node) o;
        return n.getNum() == this.num;
    }
}

class GUI extends JFrame {
    private static final int WIDTH = 960;
    private static final int HEIGHT = 960;
    private static final String TITLE = "Graph Visualisation";

    private Population initialPopulation;
    private Ordering bestOrdering;
    private NodeList nodelist;
    private int generations;

    private JPanel mainPanel;
    private JPanel bottomPanel;
    private Painter graphPainter;
    private JLabel orderingLabel;
    private JLabel currentGenerationLabel;
    private JLabel bestOrderingFitnessCost;
    private JLabel currentFitnessFunction;
    private JButton generateNextPopulationButton;
    private JButton generateFinalPopulationButton;
    private JButton swapFitnessFunction;
    private JButton resetButton;

    public GUI(Population initialPopulation, NodeList nodelist,
            int generations, ActionListener listener) {
        super(TITLE);
        this.initialPopulation = initialPopulation;
        this.nodelist = nodelist;
        setSize(WIDTH, HEIGHT);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        initMainPanel(listener);
        setVisible(true);
    }

    public void update(Ordering ordering, double fitness, int generation) {
        orderingLabel.setText("Current best ordering  " + ordering);
        currentGenerationLabel.setText("Current generation: " + generation);
        bestOrderingFitnessCost.setText("Ordering Fitness: " + fitness);
        graphPainter.setOrdering(ordering);
        repaint();
    }

    public void setFitnessFunctionText(Population.FitnessFunction fc) {
        switch(fc) {
            case TETTAMANZI:
                currentFitnessFunction.setText("Using Tettamanzi Fitness Function");
                break;

            case DEFAULT:
                currentFitnessFunction.setText("Using Default Fitness Function");
                break;
        }

    }

    private void initMainPanel(ActionListener listener) {
        mainPanel = new JPanel();
        mainPanel.setPreferredSize(new Dimension(WIDTH, HEIGHT));
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));

        orderingLabel = new JLabel("Current best ordering  "
                + initialPopulation.getBestFromGeneration());
        currentGenerationLabel = new JLabel("Current generation: 0");
        bestOrderingFitnessCost = new JLabel("Ordering Fitness: ");
        currentFitnessFunction = new JLabel("Using default fitness function");
        generateNextPopulationButton = new JButton("Next Generation");
        generateFinalPopulationButton = new JButton("Last Generation");
        swapFitnessFunction = new JButton("Swap Fitness Function");
        resetButton = new JButton("Reset");

        graphPainter = new Painter(initialPopulation.getOrdering(0), nodelist,
                initialPopulation.getChunk());

        bottomPanel = new JPanel();
        bottomPanel.setLayout(new BoxLayout(bottomPanel, BoxLayout.Y_AXIS));
        bottomPanel.setPreferredSize(new Dimension(WIDTH, 200));
        bottomPanel.add(orderingLabel);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(currentGenerationLabel);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(bestOrderingFitnessCost);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(currentFitnessFunction);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(generateNextPopulationButton);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(generateFinalPopulationButton);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(resetButton);
        bottomPanel.add(Box.createVerticalGlue());
        bottomPanel.add(swapFitnessFunction);
        bottomPanel.add(Box.createVerticalGlue());
        mainPanel.add(bottomPanel);
        mainPanel.add(graphPainter);
        add(mainPanel);

        generateNextPopulationButton.addActionListener(listener);
        generateFinalPopulationButton.addActionListener(listener);
        resetButton.addActionListener(listener);
        swapFitnessFunction.addActionListener(listener);
    }
}

class Painter extends JPanel {
    private static final int WIDTH = 960;
    private static final int HEIGHT = 700;
    private static final int RADIUS = 250;
    private static final int SHIFTX = WIDTH/2;
    private static final int SHIFTY = HEIGHT/2;

    private Ordering ordering;
    private NodeList nodelist;
    private double chunk;
    private Map<Integer, double[]> coordinates;

    public Painter(Ordering ordering, NodeList nodelist, double chunk) {
        this.ordering = ordering;
        this.nodelist = nodelist;
        this.chunk = chunk;
        coordinates = ordering.generateCoordinateMap(chunk);
        setPreferredSize(new Dimension(WIDTH, HEIGHT));
    }

    public void setOrdering(Ordering ordering) {
        this.ordering = ordering;
        coordinates = ordering.generateCoordinateMap(chunk);
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;

        g2.setStroke(new BasicStroke(2));
        g2.setColor(Color.red);
        g2.drawOval((WIDTH/2)-RADIUS, (HEIGHT/2)-RADIUS, RADIUS*2, RADIUS*2);
        g2.setColor(Color.black);

        for(int i = 0; i < ordering.length; i++) {
            Node node = nodelist.getNode(ordering.get(i));
            ArrayList<Integer> nodeConnections = node.getTails();
            double[] nodePoint = coordinates.get(node.getNum());
            //System.out.println(node.getNum() + ": x:" + nodePoint[0] + " y:" + nodePoint[1]);

            drawOrderingValue(nodePoint, i, node.getNum(), g2);
            drawNodePoint(nodePoint, g2);
            drawEdges(nodePoint, nodeConnections, g2);
        }
    }

    private void drawEdges(double[] nodePoint, ArrayList<Integer> connections,
            Graphics2D g) {
        for(Integer nodeValue : connections) {
            double[] connectionPoint = coordinates.get(nodeValue);
            g.setColor(Color.black);
            g.setStroke(new BasicStroke(3));
            g.drawLine(
                    (int)(nodePoint[0] * RADIUS) + SHIFTX,
                    (int)(-nodePoint[1] * RADIUS) + SHIFTY,
                    (int)(connectionPoint[0] * RADIUS) + SHIFTX,
                    (int)(-connectionPoint[1] * RADIUS) + SHIFTY);
        }

    }

    private void drawOrderingValue(double[] nodePoint, int index, int label,
            Graphics2D g) {
        g.setStroke(new BasicStroke(4));
        if (index < ordering.length / 2) {
            g.drawString(
                    Integer.toString(label),
                    (int)(nodePoint[0] * RADIUS) + SHIFTX + 7,
                    (int)(-nodePoint[1] * RADIUS) + SHIFTY + 7);
        } else {
            g.drawString(
                    Integer.toString(label),
                    (int)(nodePoint[0] * RADIUS) + SHIFTX - 14,
                    (int)(-nodePoint[1] * RADIUS) + SHIFTY - 14);
        }
    }

    private void drawNodePoint(double[] nodePoint, Graphics2D g) {
        g.setColor(Color.blue);
        g.fillOval(
                (int)(nodePoint[0] * RADIUS) + SHIFTX - 7,
                (int)(-nodePoint[1] * RADIUS) + SHIFTY - 7,
                14,
                14);
    }
}
