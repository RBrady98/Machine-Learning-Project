import java.util.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.IOException;
import javax.swing.JFrame;
import java.awt.Graphics;

public class main {
    public static void main(String[] args) {
        ArrayList<String> edgeStrings = getNodesFromFile("input.txt");
        NodeList nodes = generateNodeList(edgeStrings);
        System.out.println("Graph Adjacency Matrix");
        System.out.println();
        nodes.printAdjacencyMatrix();
        System.out.println();
        System.out.println();

        System.out.println("Generating random orderings");
        Population currentPopulation = new Population(15, 0, 0, 0, nodes);
        currentPopulation.generateRandomOrderings(nodes.size());
        System.out.println();
        currentPopulation.generateNextPopulation();

        Painter p = new Painter();
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
    private int populationSize;
    private int numberOfGenerations;
    private int crossoverRate;
    private int mutationRate;
    private double chunk;
    private NodeList nodelist;
    private Ordering[] orderings;

    public Population(int populationSize, int numberOfGenerations, int crossoverRate,
            int mutationRate, NodeList nodelist) {
        this.populationSize = populationSize;
        this.numberOfGenerations = numberOfGenerations;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.nodelist = nodelist;
        orderings = new Ordering[populationSize];
    }

    public void setOrderings(Ordering[] orderings) {
        this.orderings = orderings;
    }

    /**
     * generateRandomOrderings is the initial method used during generation 0 to
     * create the initial Orderings
     */
    public void generateRandomOrderings(int graphSize) {
        chunk = (2 * Math.PI) / (graphSize - 1);
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
        performSelection();

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
                    // Select indexes from number of remaining orders which still later
                    // have to be changed to actual remaining ordering indexes
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

                    orderPositionsRemaining.remove(orderRemainingIndex1);
                    orderPositionsRemaining.remove(orderRemainingIndex2);

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

        Population nextPopulation = new Population(populationSize, numberOfGenerations,
                crossoverRate, mutationRate, nodelist);
        nextPopulation.setOrderings(nextGenerationOrderings);

        System.out.println("Next Gen orderings");
        for(int i = 0; i < nextGenerationOrderings.length; i++) {
            System.out.println((i + 1)  + " => " + nextGenerationOrderings[i]);
        }

        return nextPopulation;
    }

    /**
     * calculateFitnessForAllOrderings loops through all orderings in the orderings
     * array and calculates that orderings fitness cost
     */
    private void calculateFitnessForAllOrderings() {
        for(Ordering ordering : orderings) {
            calculateOrderingFitness(ordering, chunk);
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

        Ordering[] topThird = Arrays.copyOfRange(orderings, 0,
                (int)(populationSize / 3));
        int lastThirdStartIndex = (int)((populationSize / 3) * 2);
        for(int i = lastThirdStartIndex, j = 0; i < populationSize; i++, j++) {
            orderings[i] = topThird[j];
        }
    }

    /**
     * CalculateOrderingFitness calculates the fitness for a single ordering
     * by summing the length of all edges
     */
    private void calculateOrderingFitness(Ordering ordering, double chunk) {
        double fitnessCost = 0;
        Map<Integer, double[]> coordinates = ordering.generateCoordinateMap(chunk);

        for(int i = 0; i < ordering.length; i++) {
            Node node = nodelist.getNode(ordering.get(i));
            ArrayList<Integer> nodeConnections = node.getTails();

            for(Integer nodeValue : nodeConnections) {
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
        ordering.setFitnessCost(fitnessCost);
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
        HashMap<Integer, ArrayList<Integer>> nodeTailMap
            = new HashMap<Integer, ArrayList<Integer>>();
        for(Node n : nodes) {
            nodeTailMap.put(n.getNum(), n.getTails());
        }

        //Print top bar of adjacency matrix
        System.out.printf("%3s", " |");
        nodeTailMap.forEach((k, v) -> System.out.printf("%3s", k + "|"));
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

class Painter extends JFrame {
    private static final String TITLE = "Graph Visualisation";
    private static final int WIDTH = 960;
    private static final int HEIGHT = 960;
    private Ordering ordering;
    private Map<Integer, double[]> coordinates;

    //public Painter(Ordering ordering, Map<Integer, double[]> coordinates) {
        //this.ordering = ordering;
        //this.coordinates = coordinates;
        //setTitle(TITLE);
        //setSize(WIDTH, HEIGHT);
        //setVisible(true);
        //setDefaultCloseOperation(EXIT_ON_CLOSE);
    //}

    public Painter() {
        setTitle(TITLE);
        setSize(WIDTH, HEIGHT);
        setVisible(true);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
    }

    @Override
    public void paint(Graphics g) {
        int radius = 150;

    }
}
