
import java.util.*;
import java.io.*;

public class myrecsys {

    public static String trainingData; //user id | item id | rating | timestamp
    public static String testData; //user id | item id | rating | timestamp
    public static String algorithm;
    public static HashMap<Integer, Map<Integer, Integer>> allUsers = new HashMap<Integer, Map<Integer, Integer>>();
    public static Set<Integer> userID = new HashSet<Integer>(); //keep track of all the different user id's
    public static Set<Integer> movieID = new HashSet<Integer>(); //keep track of all the different movie id's
    public static HashMap<Integer, Double> totalAvg = new HashMap<Integer, Double>(); //holds total average rating for each movie => movieID | avgRating
    public static HashMap<Integer, Set<Integer>> nonRaters = new HashMap<Integer, Set<Integer>>(); //holds all the movies with their list of users that didn't rate each => movieID | Set of all users
    public static HashMap<Integer, Map<Integer, Double>> userSimiliar = new HashMap<Integer, Map<Integer, Double>>(); //Euclidean similarities
    public static HashMap<Integer, Map<Integer, Double>> pearsonSimiliar = new HashMap<Integer, Map<Integer, Double>>(); //Pearson similarities
    public static HashMap<Integer, Set<Integer>> movieRaters = new HashMap<Integer, Set<Integer>>(); //holds movie -> all users who rated said movie
    public static HashMap<Integer, Map<Integer, Integer>> testD = new HashMap<Integer, Map<Integer, Integer>>(); //.test data set userID | <Movie, Rating>
    public static HashMap<Integer, Map<Integer, Double>> userNeighborhood = new HashMap<Integer, Map<Integer, Double>>();
    public static HashMap<Integer, Map<Integer, Double>> cosineSimiliar = new HashMap<Integer, Map<Integer, Double>>(); //cosine Similarities

    //arguments in following format: trainingData, testData, algorithm to use
    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.print("Not enough arguments entered\nNeeded arguments: trainingData file, testData file, "
                    + "algorithm to use\n");
            System.exit(-1);
        }
        trainingData = args[0];
        testData = args[1];
        algorithm = args[2]; 

        try {
            readDataBase(trainingData);
            readDataTest(testData);
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        if (algorithm.toLowerCase().equals("average")) {
            getMovieRaterPairs();
            getNonRaters();
            average();
            double error = 0;
            double size = 0;
            for (Map.Entry<Integer, Set<Integer>> e : nonRaters.entrySet()) { //movie, users
                double predicted = totalAvg.get(e.getKey());
                int m = e.getKey();
                for (int user : e.getValue()) {
                    size++;
                    if (testD.get(user) != null && testD.get(user).containsKey(m)) {
                        double actual = testD.get(user).get(m).doubleValue();
                        error += calculateRMSE(predicted, actual);
                    }

                }
            }
            error = Math.sqrt(error / size);
            writeData(error);

        } else if (algorithm.toLowerCase().equals("user-euclidean")) {
            getMovieRaterPairs();
            getNonRaters();
            userEuclidean();

        } else if (algorithm.toLowerCase().equals("user-pearson")) {
            getMovieRaterPairs();
            getNonRaters();
            userPearson(); //calculates all of the pearson coefficients between two users
        } else if (algorithm.toLowerCase().equals("item-cosine")) {
            getMovieRaterPairs();
            getNonRaters();
            itemCosine();
            neighborhood(100, cosineSimiliar);
            calculateCosineError();
        }
    }
    
//Reads in a base file of training set data and parses the information appropriately into a hashmap of all users as well as to Sets of movie ids and user ids
    public static void readDataBase(String file) throws IOException {
        HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
        Scanner reader = new Scanner(new File(file));
        //int holder = 0;
        while (reader.hasNext()) {
            String[] line = reader.nextLine().split("\\s+");
            int id = Integer.parseInt(line[0]);
            int movie = Integer.parseInt(line[1]);
            userID.add(id);
            movieID.add(movie);
            //System.out.println("movie => " + movie + " id: " + id);
            int rating = Integer.parseInt(line[2]);
            if (line.length > 3) {
                int time = Integer.parseInt(line[3]);
            }

            if (allUsers.containsKey(id)) { //first time seeing id
                map.put(movie, rating);
                allUsers.put(id, map);
            } else if (!(allUsers.containsKey(id))) {
                map = new HashMap<Integer, Integer>();
                map.put(movie, rating);
                allUsers.put(id, map);
            }
        }
    }

    //reads in a test data set and parses the information int the TestD hashmap
    public static void readDataTest(String file) throws IOException {
        //HashSet<Integer> usertemp = new HashSet<Integer>();
        HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
        Scanner reader = new Scanner(new File(file));
        //int holder = 0;
        while (reader.hasNext()) {
            String[] line = reader.nextLine().split("\\s+");
            int id = Integer.parseInt(line[0]);
            int movie = Integer.parseInt(line[1]);
            //System.out.println("movie => " + movie + " id: " + id);
            int rating = Integer.parseInt(line[2]);
            if (line.length > 3) {
                int time = Integer.parseInt(line[3]);
            }
            if (testD.containsKey(id)) { //first time seeing id
                if (rating >= 1) {
                    map.put(movie, rating);
                }
                if (id >= 0) {
                    testD.put(id, map);
                }
            } else if (!(testD.containsKey(id))) {
                map = new HashMap<Integer, Integer>();
                if (rating >= 1) {
                    map.put(movie, rating);
                }
                if (id >= 0) {
                    testD.put(id, map);
                }
            }
        }

        reader.close();
    }
    
    public static void writeData(double error) {
        System.out.println("MYRESULTS Training\t=\t<" + trainingData + ">");
        System.out.println("MYRESULTS Testing\t=\t<" + testData + ">");
        System.out.println("MYRESULTS Algorithm\t=\t<" + algorithm + ">");
        System.out.println("MYRESULTS RMSE\t\t=\t<" + error + ">");
    }
    /*----------------------------------------------------------------------------
     * Calculates the average of all the movies 
     */
    public static void average() {
        //HashMap<Integer, Double> totalAvg = new HashMap<Integer, Double>(); //takes total avg of all ratings of a given movie
        //HashSet<Integer> noUserRate = new HashSet<Integer>();
        double avg = 0;
        int totalReviewed = 0;
        int totalRate = 0;
        for (int mID : movieID) {
            for (Map.Entry<Integer, Map<Integer, Integer>> entry : allUsers.entrySet()) {
                int userid = entry.getKey();
                if (allUsers.get(userid).containsKey(mID)) {
                    totalRate += allUsers.get(userid).get(mID).intValue();
                    totalReviewed++;
                }
            }
            avg = (double) totalRate / (double) totalReviewed; //get average rating for that movie
            totalAvg.put(mID, avg); //put movie with its average into hashmap for storing
            avg = 0;
            totalReviewed = 0;
            totalRate = 0; //reset values

        }

    }
    /*--------------------------------------------------------------------------------------------------
     * Find similiarities between users and use that to predict a rating for a given user
     * Compare the predicted score with the actual score in the .test file
     */
    public static void userEuclidean() {
        getSimiliarityPairs(); //user | user2, similiarity between the two
        neighborhood(100, userSimiliar);
        calcEuclideanError();
    }
//Calculates the error based on the calculations of the similarities using the euclideanean distance metric
    public static void calcEuclideanError() {
        double error = 0;
        double total = 0, totalSim = 0, predicted = 0, actual = 0, size = 0;
        for (Map.Entry<Integer, Set<Integer>> e : nonRaters.entrySet()) { //movie, users
            int m = e.getKey();
            for (int user : e.getValue()) { //cycle through users who did not rate the movie
                size++;
                for (Map.Entry<Integer, Map<Integer, Double>> entry : userNeighborhood.entrySet()) { //cycle through users
                    int utemp = entry.getKey();
                    if (movieRaters.get(m).contains(utemp)) { //if movie was rated by user
                        if (userNeighborhood.containsKey(user)) {
                            if (userNeighborhood.get(user).containsKey(utemp)) { //if the two users share a similarity in the neighborhood
                                double sim = userNeighborhood.get(user).get(utemp).doubleValue();
                                double temprate = allUsers.get(utemp).get(m).doubleValue();
                                total += (temprate * sim);
                                if (testD.get(user) != null && testD.get(user).containsKey(m)) {
                                    totalSim += sim;
                                    actual = testD.get(user).get(m).doubleValue(); //get actual score                          
                                }
                            }

                        }

                    }
                }
                if (totalSim != 0) {
                    predicted = total / totalSim;
                } else {
                    predicted = 0;
                }
                totalSim = 0;
                total = 0;
                error += calculateRMSE(predicted, actual);
            }
        }
        error = Math.sqrt(error / size);

        writeData(error);
    }

    /*----------------------------------------------------------------------------------------------------------
     * Gets all the movies with their rated users with their ratings
     */
    public static void getMovieRaterPairs() {
        for (int mID : movieID) {
            Set<Integer> temp = new HashSet<Integer>();
            for (Map.Entry<Integer, Map<Integer, Integer>> entry : allUsers.entrySet()) {
                int userid = entry.getKey();
                if (allUsers.get(userid).containsKey(mID)) {
                    temp.add(userid);
                }
            }
            movieRaters.put(mID, temp);
            //temp.clear();
        }
    }
//Gets all the users that didn't rate a particular movie
    public static void getNonRaters() {
        for (int mID : movieID) {
            Set<Integer> temp = new HashSet<Integer>();
            for (Map.Entry<Integer, Map<Integer, Integer>> entry : allUsers.entrySet()) {
                int id = entry.getKey();
                if (!(movieRaters.get(mID).contains(id))) {
                    temp.add(id);
                }
            }
            nonRaters.put(mID, temp);
        }
    }
    /*------------------------------------------------------------------------------------------------------------
     * Calculates all the similiarity pairs between the users based on the getDistance and getSimiliarity methods
     * Similiarity = 1 / (1+distance(i,j))
     * Adds to the global HashMap containing user1 as key and user2 and their similiarity as another map
     */
    public static void getSimiliarityPairs() {
        for (Map.Entry<Integer, Map<Integer, Integer>> entry : allUsers.entrySet()) {
            HashMap<Integer, Double> maptemp = new HashMap<Integer, Double>();
            int user1 = entry.getKey();
            for (Map.Entry<Integer, Map<Integer, Integer>> e : allUsers.entrySet()) {
                int user2 = e.getKey();
                if (!(user1 == user2)) {
                    double s = calcSimiliarity(allUsers.get(user1), allUsers.get(user2));
                    maptemp.put(user2, s);
                    userSimiliar.put(user1, maptemp);
                }
            }
        }
    }
    /*-----------------------------------------------------------------------------------
     * Calculates similiarity between two users
     * Formula: 1/(1+distance(user1, user2)
     * Takes in two users and returns a double representing their similiarity
     */
    public static double calcSimiliarity(Map<Integer, Integer> user1, Map<Integer, Integer> user2) {
        double similiarity = 1 / (1 + calcDistance(user1, user2));
        if (1 + calcDistance(user1, user2) == 0) {
            return 0;
        }
        return similiarity;
    }

    /*-----------------------------------------------------------------------------------
     * Calculates the distance between two users based on differences in ratings of 
     * movies both have rated
     * Takes in two users
     * Returns a double corresponding to the distance between the two users
     */
    public static double calcDistance(Map<Integer, Integer> user1, Map<Integer, Integer> user2) {
        int N = 0;
        double totalDifference = 0;
        HashMap<Integer, Integer> u1 = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> u2 = new HashMap<Integer, Integer>();
        for (Map.Entry<Integer, Integer> entry : user1.entrySet()) {
            u1.put(entry.getKey(), entry.getValue());
        }
        for (Map.Entry<Integer, Integer> entry : user2.entrySet()) {
            u2.put(entry.getKey(), entry.getValue());
        }
        for (int mID : movieID) {
            //System.out.println(u1.containsKey(mID) && u2.containsKey(mID));
            if (u1.containsKey(mID) && u2.containsKey(mID)) {
                N++;
                int r1 = u1.get(mID);
                int r2 = u2.get(mID);
                totalDifference += ((r1 - r2) * (r1 - r2));
            }
        }
        //System.out.println(totalDifference + "<- difference  N-> " + N);
        return Math.sqrt(totalDifference);
    }

    //returns pearson coefficient between two users based on movie ratings
    public static double getPearson(double[] x, double[] y) {
        double x_bar = 0;
        double y_bar = 0;
        for (int i = 0; i < x.length; i++) {
            x_bar += x[i];
            y_bar += y[i];
        }
        x_bar /= x.length;
        y_bar /= y.length;
        double totalTop = 0;
        double totalBottomx = 0, totalBottomy = 0;
        for (int i = 0; i < x.length; i++) {
            totalTop += (x[i] - x_bar) * (y[i] - y_bar);
        }
        for (int i = 0; i < x.length; i++) {
            totalBottomx += (Math.pow((double) (x[i] - x_bar), 2));
            totalBottomy += (Math.pow((double) (y[i] - y_bar), 2));
        }
        if (Math.sqrt(totalBottomx) * Math.sqrt(totalBottomy) == 0) {
            return 0;
        }
        return totalTop / (Math.sqrt(totalBottomx) * Math.sqrt(totalBottomy));
    }

    //Selects k users that have the highest similiarity with the active user (neighborhood), set to 100 for all algorithms
    //applied to the euclidean, pearson, and cosine algorithms
    public static void neighborhood(int k, HashMap<Integer, Map<Integer, Double>> map) {
        int count = 0;
        ArrayList<Double> list = new ArrayList<Double>();
        for (Map.Entry<Integer, Map<Integer, Integer>> entry : allUsers.entrySet()) {
            int user1 = entry.getKey();
            if (!map.containsKey(user1)) {
                continue;
            }
            HashMap<Integer, Double> t = new HashMap<Integer, Double>();
            list = new ArrayList<Double>();
            ArrayList<Integer> u = new ArrayList<Integer>();
            Map.Entry<Integer, Double> min = null;
            double minimum = 0;
            for (Map.Entry<Integer, Double> e : map.get(user1).entrySet()) {
                if (min == null || min.getValue() > e.getValue()) {
                    min = e;
                    minimum = e.getValue();
                }
                if (count < k) { //if we are looking at the first k users just add to list
                    list.add(e.getValue());
                    u.add(e.getKey());
                    count++;
                } else { //we have more than k users now we need to see if this similiarity is higher than our lowest
                    if (e.getValue() > minimum || e.getValue() > Collections.min(list)) {
                        int index = list.indexOf(Collections.min(list));
                        list.remove(Collections.min(list));
                        u.remove(index);
                        list.add(e.getValue());
                        u.add(e.getKey());
                    }
                }
            }
            for (int i = 0; i < list.size(); i++) {
                t.put(u.get(i), list.get(i));
            }
            count = 0;
            userNeighborhood.put(user1, t);
        }

    }

    /*---------------------------------------------------------------------------------
     * r = (n(sum of x*y)-(sum of x)*(sum of y))/(sqrt((n*sum of (x^2)-(sum of x)^2)*(n*sum of (y^2)-(sum of y)^2))
     * General case:
     * 1. Collect ratings from users
     * 2. assign a weight to all users with respect to similarity with the active user
     * 3. Select k users that have the highest similarity with the active user (neighborhood) (=100)
     * 4. Compute a prediction score from a weighted combination of the selected neighborhood ratings
     */
    public static void userPearson() {
        int n = 0;
        double r = 0;
        HashMap<Integer, Double> temp = new HashMap<Integer, Double>();
        for (int user1 : userID) {
            temp = new HashMap<Integer, Double>();
            for (int user2 : userID) {
                if (user1 != user2) {
                    ArrayList<Double> x = new ArrayList<Double>();
                    ArrayList<Double> y = new ArrayList<Double>();
                    for (Map.Entry<Integer, Integer> e1 : allUsers.get(user1).entrySet()) {
                        int movie1 = e1.getKey();
                        double mRating1 = e1.getValue();
                        if (allUsers.get(user2) != null && allUsers.get(user2).containsKey(movie1)) {
                            n++; //total of same reviewed movies incremented
                            x.add(normalize(mRating1)); //normalize rankings from [1,5] to [-1,1] for calculation
                            y.add(normalize(allUsers.get(user2).get(movie1)));
                        }
                    }
                    double[] x1 = new double[x.size()];
                    double[] y1 = new double[y.size()];
                    for (int i = 0; i < x.size(); i++) {
                        x1[i] = x.get(i);
                        y1[i] = y.get(i);
                    }
                    r = getPearson(x1, y1);
                    temp.put(user2, r);
                    pearsonSimiliar.put(user1, temp);
                }
            }
        }
        neighborhood(100, pearsonSimiliar);
        calculatePearsonError();
    }

    //Calculates the error associated with the given similarities generated and the difference between the weighted prediction vs the actual rating
    public static void calculatePearsonError() {
        double total = 0, totalSim = 0, predicted = 0, actual = 0, size = 0;
        double error = 0;
        for (Map.Entry<Integer, Set<Integer>> e : nonRaters.entrySet()) { //movie, users
            int m = e.getKey();
            for (int user : e.getValue()) { //cycle through users who did not rate the movie
                size++;
                for (Map.Entry<Integer, Map<Integer, Double>> entry : userNeighborhood.entrySet()) { //cycle through users
                    int utemp = entry.getKey();
                    if (movieRaters.get(m).contains(utemp)) {
                        if (userNeighborhood.containsKey(user)) {
                            if (userNeighborhood.get(user).containsKey(utemp)) {
                                double sim = (userNeighborhood.get(user).get(utemp).doubleValue());
                                double temprate = normalize(allUsers.get(utemp).get(m).doubleValue());
                                total += temprate * sim;
                                if (testD.get(user) != null && testD.get(user).containsKey(m)) {
                                    totalSim += Math.abs(sim);
                                    actual = testD.get(user).get(m).doubleValue(); //get actual score                          
                                }
                            }

                        }

                    }
                }
                if (totalSim != 0) {
                    predicted = total / totalSim;
                } else {
                    predicted = 0;
                }
                totalSim = 0;
                total = 0;
                error += calculateRMSE(denormalize(predicted), actual);
            }
        }
        error = Math.sqrt(error / size);
        writeData(error);
    }
    
    /*-------------------------------------------------------------------------------------------
     * similarity = cos(O) = (A dot B)/(||A||*||B||)  == (sum of Ai x Bi)/(sqrt(sum of (Ai)^2) x sqrt(sum of (Bi)^2))
     */
    public static void itemCosine() {
        double c = 0;
        HashMap<Integer, Double> temp = new HashMap<Integer, Double>();
        for (int m1 : movieID) {
            temp = new HashMap<Integer, Double>();
            for (int m2 : movieID) {
                if (m1 != m2) {
                    ArrayList<Double> x = new ArrayList<Double>();
                    ArrayList<Double> y = new ArrayList<Double>();
                    for (Integer e1 : movieRaters.get(m1)) {
                        int user1 = e1;
                        double mRating1 = allUsers.get(user1).get(m1).doubleValue();
                        if (allUsers.get(user1) != null && allUsers.get(user1).containsKey(m2)) {
                            x.add((mRating1));
                            y.add((allUsers.get(user1).get(m2).doubleValue()));
                        }
                    }
                    double[] x1 = new double[x.size()];
                    double[] y1 = new double[y.size()];
                    for (int i = 0; i < x.size(); i++) {
                        x1[i] = x.get(i);
                        y1[i] = y.get(i);
                    }
                    c = calculateCosine(x1, y1);
                    temp.put(m2, c);
                    cosineSimiliar.put(m1, temp);
                }
            }
        }
    }
    //By user add up ratings of user1 * sim of movie1 to movie2
    public static void calculateCosineError() {
        double total = 0, totalSim = 0, predicted = 0, actual = 0, size = 0;
        double error = 0;
        for (Map.Entry<Integer, Set<Integer>> e : nonRaters.entrySet()) { //movie, users
            int m1 = e.getKey();
            for (int user : e.getValue()) { //cycle through users who did not rate the movie
                size++;
                for (Integer mID2 : movieID) { // m1 | m2 | similarity
                    if (m1 == mID2) {
                        continue;
                    }
                    if (movieRaters.get(mID2).contains(user)) {
                        if (userNeighborhood.containsKey(m1)) {
                            if (userNeighborhood.get(m1).containsKey(mID2)) {
                                double sim = userNeighborhood.get(m1).get(mID2).doubleValue();
                                double temprate = (allUsers.get(user).get(mID2).doubleValue());
                                total += (temprate * sim);
                                if (testD.get(user) != null && testD.get(user).containsKey(m1)) {
                                    totalSim += Math.abs(sim);
                                    actual = testD.get(user).get(m1).doubleValue(); //get actual score                          
                                }
                            }

                        }
                    }
                }
                if (totalSim != 0) {
                    predicted = total / totalSim;
                } else {
                    predicted = 0;
                }
                totalSim = 0;
                total = 0;
                error += calculateRMSE((predicted), actual);
            }
        }
        error = Math.sqrt(error / size);
        writeData(error);
    }

    //calculates the cosine similarity based of of two arrays symbolizing the vectors
    public static double calculateCosine(double[] x, double[] y) {
        //= A.B/ |A|*|B|
        double totalTop = 0, totalX = 0, totalY = 0;
        double totalBottom = 0;
        for (int i = 0; i < x.length; i++) {
            totalTop += x[i] * y[i];
            totalX += x[i] * x[i];
            totalY += y[i] * y[i];
        }
        totalBottom = Math.sqrt(totalX) * Math.sqrt(totalY);
        if (totalBottom == 0) {
            return 0;
        }
        return (totalTop / totalBottom);
    }

    //Returns the normalized ranking of [1,5] to [-1,1]
    public static double normalize(double rank) {
        return ((.5 * (double) rank) - 1.5);
    }
    
    //Returns denormalized ranking of [-1, 1] to [1,5]
    public static double denormalize(double nRank) {
        return (2 * nRank + 3);
    }

    /*------------------------------------------------------------------------------------------
     * RMSE = sqrt((sum of (pi-qi)^2)/N)
     */
    public static double calculateRMSE(double predicted, double actual) {
        double rmse = 0;
        double diff = 0;
        diff = predicted - actual;
        rmse = diff * diff;
        return rmse;
    }

}
