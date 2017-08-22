import com.google.gson.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class ProcessTwitterData {

    private static final String TWITTER="EEE MMM dd HH:mm:ss ZZZZZ yyyy";
    final static long MONTH = 2629743;
    final static long WEEK = 604800;
    final static long DAY = 86400;

    private long timeAggregation;

    String path;
    HashMap<Integer, HashMap<Long, Set<String>>> userHashtags;
    List<String> hashlist;

    String[] hashtags;

    public ProcessTwitterData(String path, long timeAggregation) {
        this.timeAggregation = timeAggregation;
        if (path == null) this.path = "./userTweets/";

        userHashtags = new HashMap<>();
        File folder = new File(path);
        File[] listOfFiles = folder.listFiles();

        if (listOfFiles.length == 0) System.err.println("NO FILES FOUND");
        else {
            Arrays.stream(listOfFiles)
                .filter(thing -> thing.getName().endsWith(".txt"))
                .forEach(file -> processFile(file));
        }
        System.out.println(userHashtags);
        mapHashtags("HashtagMapping.txt");
        processData();
    }

    private void processFile(File file) {
        try (Stream<String> stream = Files.lines(file.toPath())) {

//            stream.forEach(System.out::println);
            stream.forEach(this::processJSON);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void processJSON(String json) {

        JsonArray array = new JsonParser().parse(json).getAsJsonArray();
        for (int i = 0; i<array.size(); i++){
            JsonObject tweet = array.get(i).getAsJsonObject();
            if (tweet.get("text").toString().chars().filter(ch -> ch == '#').count() > 1) {
                Integer userId = tweet.getAsJsonObject("user").get("id").getAsInt();
                Long epochDate = parseDate(tweet.get("created_at").getAsString());
//                int date = aggregate(epochDate);
                Set<String> hashtags = getHashTags(tweet.get("text").getAsString());

                    if (userHashtags.containsKey(userId)) {
                        if (userHashtags.get(userId).containsKey(epochDate)) {
                            userHashtags.get(userId).get(epochDate).addAll(hashtags);
                        } else {
                            userHashtags.get(userId).put(epochDate, hashtags);
                        }
                    } else {
                        HashMap<Long, Set<String>> dateHashtags = new HashMap<>();
                        dateHashtags.put(epochDate, hashtags);
                        userHashtags.put(userId, dateHashtags);
                    }
                }
            }
    }

    private Set<String> getHashTags(String text) {
        Set<String> result = new HashSet<>();
        Pattern pattern = Pattern.compile("#(\\S+)");
        Matcher mat = pattern.matcher(text);
        while(mat.find()) {
            result.add(mat.group(1));
        }
        return result;
    }

    private long parseDate(String date) {

        SimpleDateFormat sf = new SimpleDateFormat(TWITTER);
        sf.setLenient(true);
        try {
            return sf.parse(date).getTime() / 1000;
        } catch (ParseException e) {
            System.err.println("Failed to parse date");
            e.printStackTrace();
        }
        return 0;
    }

    private int aggregate(long date) {
        return Math.toIntExact(date / this.timeAggregation);
    }

    private void mapHashtags(String filename){

        // First assign a number to each hashtag
        Set<String> allHashtags = new HashSet<>();
        userHashtags.values().parallelStream().map(map -> map.values()).forEach(set -> set.forEach(tag -> allHashtags.addAll(tag)));
        hashtags = allHashtags.toArray(new String[0]);
        System.out.println("Total hashtags: "+hashtags.length);

        try {
            PrintWriter writer = new PrintWriter(filename);
            writer.println("Hashtag\tID");
            for(int i = 0; i<hashtags.length; i++) {
                writer.println(hashtags[i] + "\t" + i);
            }
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private Set<Edge> allHashtagPairs(Set<String> hashtags) {
        Set<Edge> result = new HashSet<>();
        for(String s1 : hashtags) {
            for(String s2 : hashtags) {
                if(!s1.equals(s2)) {
                    int node1 = hashlist.indexOf(s1);
                    int node2 = hashlist.indexOf(s2);
                    result.add(new Edge(node1, node2));
                }
            }
        }
        return result;
    }

    private void processData(){
        if(userHashtags == null) {
            System.out.println("No DATA");
            return;
        }

        this.hashlist = Arrays.asList(this.hashtags);

        TreeMap<Integer, HashMap<Edge, Integer>> twitterData = new TreeMap<>();

        for(Map.Entry<Integer, HashMap<Long, Set<String>>> userTweets : userHashtags.entrySet()) {
            for(Map.Entry<Long, Set<String>> timeTags : userTweets.getValue().entrySet()) {
                int userid = userTweets.getKey(); // not needed?
                int day = aggregate(timeTags.getKey());
                for (Edge e : allHashtagPairs(timeTags.getValue())) {
                    if(twitterData.containsKey(day)){
                        twitterData.get(day).compute(e, (k, v) -> v == null ? 1 : v + 1);
                    } else {
                        HashMap<Edge, Integer> hashtagWeights = new HashMap<>();
                        hashtagWeights.put(e, 1);
                        twitterData.put(day, hashtagWeights);
                    }
                }
            }
        }

        System.out.println(twitterData);
        printData(twitterData);
    }

    private void printData(TreeMap<Integer, HashMap<Edge, Integer>> twitterData) {
        try {
            PrintWriter writer = new PrintWriter("twitter_data.txt");
            writer.println("Node1,Node2,Time,Weight");
            for(int time : twitterData.keySet()) {
                for(Edge e : twitterData.get(time).keySet()) {
                    writer.println(e.getNode1()+","+e.getNode2()+","+time+","+twitterData.get(time).get(e));
                }
            }
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new ProcessTwitterData("/users/gaurav/userTweets/", ProcessTwitterData.DAY);
    }
}
